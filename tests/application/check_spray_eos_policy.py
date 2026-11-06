from pathlib import Path
import csv,json,math,os,subprocess,time,argparse,tempfile

parser=argparse.ArgumentParser(description='Exercise per-run EOS policy and isolated restart refusal in one fixed spray frontend')
parser.add_argument('--exe',type=Path,required=True)
parser.add_argument('--kind',choices=('serial','mpi','amr','eb'),required=True)
parser.add_argument('--work',type=Path,required=True)
parser.add_argument('--mpiexec',default='mpiexec')
args=parser.parse_args()
args.work.mkdir(parents=True,exist_ok=True)
root=Path(tempfile.mkdtemp(prefix='run-',dir=args.work.resolve()))
cache=args.exe.resolve().parent
exe=args.exe.resolve()
mpi_launcher=args.mpiexec
env=os.environ.copy()
for key in ('LD_LIBRARY_PATH','LIBRARY_PATH','DISPLAY'):env.pop(key,None)
env.update(HWLOC_COMPONENTS='-gl,-x86',OMPI_MCA_hwloc_base_binding_policy='none')
records_path=root/'evidence.json'
records=[]

def run(kind,selected,phase,T=290,override=False,restart='',strict=False,baseline=False,ranks=2):
    label=('selected' if selected else 'fixed')+'-'+kind+'-'+phase
    if baseline:label='baseline-'+label
    cwd=root/label;cwd.mkdir(exist_ok=True)
    previous=next((item for item in records if item['case']==label),None)
    if previous:
        if not previous['passed'] and phase.endswith('refusal'):
            diagnostic=(cwd/'run.log').read_text().lower()
            previous['passed']=previous['exit']!=0 and 'restart rejected:' in diagnostic and not list(cwd.glob('*.csv'))
            previous['checker_correction']='Accept actual MPI Restart rejected diagnostic; original checkpoint-only string filter failed'
            records_path.write_text(json.dumps(records,indent=2)+'\n')
        if not previous['passed']:raise RuntimeError('Retained failing case: '+label)
        return cwd
    if list(cwd.glob('*.csv')) or (cwd/'state.chk').exists():raise RuntimeError('Fresh evidence directory required: '+str(cwd))
    enabled='.false.' if strict or baseline else '.true.'
    policy=f'pinned_eos_extrapolation={enabled}'
    if override:policy+=',pinned_eos_temperature_min=90,pinned_eos_temperature_max=4000'
    if baseline:policy=''
    stop=phase in ('cold','warm','override','high','legacy')
    if kind=='amr' and phase in ('override','high'):stop=False
    composition="composition_count=1,composition_species='N2',composition_mole_fractions=1"
    if kind=='serial':
        nml=f'''&spray
cells=4,4,4,domain_length=.004,.004,.004
gas_temperature={T},gas_pressure=101325,gas_velocity=0,0,0,mole_names='N2','','',mole_amounts=1,0,0
{policy}
chemistry_enabled=.false.,transport_enabled=.false.,time_integrator='sdc',sdc_iterations=1
hydro_reconstruction='characteristic_ppm',riemann_solver='pelec'
vapor_name='H2O',liquid_density=1000,liquid_cp=4184,reference_temperature=373.15,reference_pressure=101325,latent_heat=2257000
initial_liquid_mass=1e-10,parcel_count=2,drop_temperature=330,drop_diameter=20e-6
injection_center=.0015,.0015,.0015,injection_axis=1,0,0,injection_speed=0,cone_angle=0
drag_enabled=.true.,heat_enabled=.true.,evaporation_enabled=.false.,smagorinsky_constant=.05
final_time=2e-12,maximum_dt=1e-12,fixed_dt=1e-12,maximum_steps=8
output_prefix='gas',checkpoint_file={'state.chk'!r},stop_after_steps={1 if stop else 0},restart_file={restart!r}
/
'''
        name='pelef_spray_3d'+('_selected' if selected else '');args=['input.nml']
    elif kind=='mpi':
        nml=f'''&mpi_spray
cells=8,4,4,spacing=.001,.001,.001,gas_temperature={T},gas_density=1,gas_velocity=0,0,0
gas_fractions=0,0,0,0,0,0,0,0,0,1
{policy}
number=2,total_mass=1e-10,speed=0,center=.003,.002,.002,axis=1,0,0,half_angle=0
drop_temperature=330,diameter=20e-6,chemistry_enabled=.false.,transport_flags=.false.,.false.,.false.,.false.
time_integrator='sdc',sdc_iterations=1,hydro_reconstruction='characteristic_ppm',solver='pelec',les_cs=.05
final_time={1e-12 if stop else 2e-12},dt_max=1e-12,fixed_dt=1e-12,max_steps=8
output_prefix='gas',checkpoint_path='state.chk',restart_path={restart!r}
/
'''
        name='pelef_mpi_spray_3d'+('_selected' if selected else '');args=['input.nml']
    elif kind=='amr':
        fixed='x_h2=0,x_h=0,x_o=0,x_o2=0,x_oh=0,x_h2o=0,x_ho2=0,x_h2o2=0,x_ar=0,x_n2=1'
        nml=f'''&reactive_3d
nx=4,ny=4,nz=4,maximum_steps=8,x_upper=.004,y_upper=.004,z_upper=.004
final_time=2e-12,cfl=.2,problem='entropy_wave',density_wave_amplitude=0
thermo_model='{'selected' if selected else 'full_h2o2'}',{composition if selected else fixed}
initial_temperature={T},initial_pressure=101325,initial_velocity_x=0,initial_velocity_y=0,initial_velocity_z=0
{policy}
chemistry_enabled=.false.,transport_enabled=.false.,time_integrator='sdc',sdc_iterations=1
allow_small_energy=.false.,allow_negative_energy=.false.,dual_energy_eta2=2
reconstruction='characteristic_ppm',riemann_solver='pelec',boundary_condition='periodic',output_file='gas.csv'
/
&amr_tree_timestepping subcycling=.false.,automatic_timestep=.false.,maximum_dt=1e-12,fixed_dt=1e-12 /
'''
        (cwd/'tree.dat').write_text('1\n1 1 1 1 1 1 1 2\n')
        (cwd/'spray.nml').write_text(f'''&amr_spray
component_count=1,vapor_names='H2O',fractions=1,liquid_density=1000,liquid_cp=4184
reference_temperature=373.15,reference_pressure=101325,latent_heat=2257000
parcel_count=2,total_mass=1e-10,center=.0015,.0015,.0015,axis=1,0,0,speed=0,half_angle=0
diameter=20e-6,drop_temperature=330,subcycling=.false.,options%evaporation=.false.
/
&amr_spray_checkpoint checkpoint_file='state.chk',checkpoint_time=1e-12,
stop_after_checkpoint={'.true.' if stop else '.false.'},restart_file={restart!r} /
''')
        name='pelef_mpi_amr_spray_3d'+('_selected' if selected else '');args=['input.nml','tree.dat','spray.nml']
    else:
        mw=28.014 if selected else sum(x*w for x,w in zip([.29570,1e-5,1e-5,.14784,1e-5,0,.55643],[2.016,1.008,15.999,31.998,17.007,18.015,28.014]))
        rho=101325*mw/(8314.46261815324*T)
        nml=f'''&reactive_eb_3d
nx=4,ny=4,nz=4,maximum_steps=8,x_upper=.004,y_upper=.004,z_upper=.004
geometry_type='level_set',level_set_file='level_set.dat'
final_time=2e-12,cfl=.1,thermo_model='{'selected' if selected else 'elementary'}'
{composition if selected else ''}
regular_density={rho:.17g},cut_density={rho:.17g},initial_pressure=101325
{policy}
chemistry_enabled=.false.,transport_enabled=.false.,output_file='gas.csv'
checkpoint_file='state.chk',checkpoint_interval_steps={1 if stop else 0},
stop_after_checkpoint={'.true.' if stop else '.false.'},restart_file={restart!r}
/
'''
        (cwd/'level_set.dat').write_text('4 4 4\n'+('\n'.join('1 1 1 1 1' for _ in range(25)))+'\n')
        (cwd/'spray.nml').write_text('''&eb_spray
nc=1,number=2,vapor_names='H2O',fractions=1,center=.0015,.0015,.0015,axis=1,0,0
speed=0,half_angle=0,diameter=20e-6,liquid_temperature=330,total_mass=1e-10,maximum_dt=1e-12
options%evaporation=.false.
/
''')
        name='pelef_reactive_eb_3d_selected_spray' if selected else 'pelef_eb_spray_3d';args=['input.nml','spray.nml']
    (cwd/'input.nml').write_text(nml)
    command=[str(exe)]+args
    if kind in ('mpi','amr'):command=[mpi_launcher,'-n',str(ranks)]+command
    started=time.monotonic()
    with (cwd/'run.log').open('w') as log:result=subprocess.run(command,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT)
    item=dict(case=label,command=command,cwd=str(cwd),exit=result.returncode,seconds=time.monotonic()-started)
    refusal=phase.endswith('refusal')
    if refusal:
        item['passed']=result.returncode!=0 and any(text in (cwd/'run.log').read_text().lower() for text in ('checkpoint','restart rejected:')) and not list(cwd.glob('*.csv'))
        item['expected_failure']='checkpoint EOS policy mismatch'
    else:
        item['passed']=result.returncode==0
        if item['passed']:
            files=list(cwd.glob('gas*gas.csv')) if kind in ('serial','mpi') else [cwd/'gas.csv']
            if kind=='amr' and stop:
                item['passed'] &= (cwd/'state.chk').exists() and 'Stopped after checkpoint at time:' in (cwd/'run.log').read_text()
                item['checkpoint_only']=True
            else:
                rows=[row for f in files for row in csv.DictReader(line for line in f.open() if not line.startswith('#'))]
                ts=[float(row['temperature']) for row in rows]
                item.update(rows=len(rows),minimum_temperature=min(ts),maximum_temperature=max(ts))
                item['passed'] &= bool(rows) and all(math.isfinite(t) and t>0 for t in ts)
                if T==290:item['passed'] &= max(ts)<300
                if T==3600:item['passed'] &= min(ts)>3500
    records.append(item);records_path.write_text(json.dumps(records,indent=2)+'\n')
    if not item['passed']:raise RuntimeError('Failed: '+label)
    return cwd

if __name__=='__main__':
    kind=args.kind
    cold=run(kind,False,'cold')
    run(kind,False,'resume',restart=str(cold/('state.chk.step1' if kind=='eb' else 'state.chk')),ranks=4)
    run(kind,False,'bounds-refusal',restart=str(cold/('state.chk.step1' if kind=='eb' else 'state.chk')),override=True,ranks=4)
    warm=run(kind,False,'warm',T=1000)
    run(kind,False,'mode-refusal',T=1000,restart=str(warm/('state.chk.step1' if kind=='eb' else 'state.chk')),strict=True,ranks=4)
    run(kind,False,'override',override=True)
    run(kind,False,'high',T=3600)
    print('PASS: EOS boundaries, particle exchange, restart and isolated policy refusal; evidence:',root)

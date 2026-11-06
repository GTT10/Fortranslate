"""A rounded final clock completes; a physically incomplete prefix still fails."""
from pathlib import Path
import argparse,csv,json,subprocess,hashlib,sys,tempfile
p=argparse.ArgumentParser()
for key in ['serial','mpi','mpiexec','work']:p.add_argument('--'+key,required=True)
a=p.parse_args();w=Path(a.work).resolve();w.mkdir(parents=True,exist_ok=True);w=Path(tempfile.mkdtemp(prefix="run_",dir=w))
base="""&eb_tree
 root_shape=4,4,4,child_count=1,parent(1)=1,ratio(1)=2,
 lower(:,1)=1,1,1,upper(:,1)=2,2,2,
 geometry_type='plane',plane_axis='z',plane_position=.31,
 dual_energy_enabled=.true.,dual_energy_eta2=2,allow_negative_energy=.false.,
 time_integrator='sdc',mol_mode=.true.,subcycling=.true.,sdc_iterations=2,
 hydro_method='plm_ctu',riemann_solver='pelec',redistribution='none',boundary_faces=6*'outflow',
 regular_density=.31,cut_density=.31,initial_pressure=250000,
 time_step=1.25e-8,final_time=2.e-7,
 chemistry_enabled=.false.,transport_enabled=.false.,
 composition_count=3,composition_species='H2','O2','N2',composition_mole_fractions=2,1,3.76,
 output_interval_steps=1,
"""
runs=[]
def run(name,steps,ranks=0,restart=None,fail=False):
 d=w/name;d.mkdir();extra=f' maximum_steps={steps},output_file="{d}/result.csv",\n'
 if restart:extra+=f' restart_file="{restart}",\n'
 else:extra+=f' checkpoint_interval_steps=16,checkpoint_file="{d}/state.chk",\n'
 (d/'input.nml').write_text(base+extra+'/\n');cmd=([str(Path(a.serial).resolve())] if not ranks else [a.mpiexec,'-n',str(ranks),str(Path(a.mpi).resolve())])+[str(d/'input.nml')]
 r=subprocess.run(cmd,cwd=d,capture_output=True,text=True,timeout=90);(d/'run.log').write_text(r.stdout+r.stderr);runs.append({'name':name,'command':cmd,'returncode':r.returncode})
 (w/'record.json').write_text(json.dumps({'status':'running','runs':runs},indent=2)+'\n')
 if fail:assert r.returncode!=0 and 'Maximum steps reached before final time' in r.stdout+r.stderr;return d
 assert r.returncode==0,(name,r.stdout+r.stderr)
 with (d/'result.csv.step16.csv.time.csv').open() as f:clock=list(csv.DictReader(f))
 assert len(clock)==1 and int(clock[0]['step'])==16
 assert abs(float(clock[0]['time_s'])-2.e-7)<=64*sys.float_info.epsilon*2.e-7
 assert not (d/'result.csv.step17.csv').exists()
 return d
serial=run('serial',16);mpi=run('mpi2',16,2)
assert (serial/'result.csv').read_bytes()==(mpi/'result.csv').read_bytes()
checkpoint=serial/'state.chk';digest=hashlib.sha256(checkpoint.read_bytes()).hexdigest();restarted=run('mpi3_completed_restart',16,3,checkpoint)
assert (serial/'result.csv').read_bytes()==(restarted/'result.csv').read_bytes()
assert hashlib.sha256(checkpoint.read_bytes()).hexdigest()==digest
run('incomplete_serial',15,fail=True);run('incomplete_mpi2',15,2,fail=True)
(w/'record.json').write_text(json.dumps({'status':'passed','cases':5,'runs':runs,'checkpoint_sha256':digest},indent=2)+'\n')

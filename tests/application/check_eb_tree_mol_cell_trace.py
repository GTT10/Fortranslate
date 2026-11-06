"""Check opt-in serial MOL diagnostics without changing ordinary trajectories."""
from pathlib import Path
import argparse,csv,json,math,os,subprocess,tempfile,collections
p=argparse.ArgumentParser()
for name in ['serial','mpi','mpiexec','work']:p.add_argument('--'+name,required=True)
a=p.parse_args();Path(a.work).mkdir(parents=True,exist_ok=True)
w=Path(tempfile.mkdtemp(prefix='mol_cell_trace_',dir=Path(a.work).resolve()))
base="""&eb_tree
 root_shape=4,4,4,child_count=1,parent(1)=1,ratio(1)=2,
 lower(:,1)=1,1,1,upper(:,1)=2,2,2,
 geometry_type='plane',plane_axis='z',plane_position=.31,
 dual_energy_enabled=.true.,dual_energy_eta2=2,allow_negative_energy=.false.,
 time_integrator='sdc',mol_mode=.true.,subcycling=.true.,sdc_iterations=2,
 hydro_method='plm_ctu',riemann_solver='pelec',redistribution='none',boundary_faces=6*'outflow',
 regular_density=.31,cut_density=.35,initial_pressure=250000,
 time_step=1.e-9,final_time=1.e-9,maximum_steps=1,
 chemistry_enabled=.false.,transport_enabled=.false.,
 composition_count=3,composition_species='H2','O2','N2',composition_mole_fractions=2,1,3.76,
 output_interval_steps=1,
"""
trace=' mol_trace_count=3,mol_trace_nodes=2,2,1,mol_trace_cells(:,1)=2,2,3,\n'+' mol_trace_cells(:,2)=2,2,4,mol_trace_cells(:,3)=3,3,3,\n'
runs=[]
def run(name,extra='',ranks=0,error=None,reacting=False):
 d=w/name;d.mkdir();inp=d/'input.nml';inp.write_text((base.replace('chemistry_enabled=.false.','chemistry_enabled=.true.') if reacting else base)+extra+f" output_file='{d}/result.csv',\n/\n")
 cmd=([str(Path(a.serial).resolve())] if not ranks else [a.mpiexec,'-n',str(ranks),str(Path(a.mpi).resolve())])+[str(inp)]
 with (d/'run.log').open('w') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,cwd=d,timeout=60)
 runs.append({'name':name,'command':cmd,'returncode':r.returncode})
 if error:assert r.returncode!=0 and error in (d/'run.log').read_text(),(name,(d/'run.log').read_text())
 else:assert r.returncode==0,(name,(d/'run.log').read_text())
 return d
baseline=run('disabled');enabled=run('enabled',trace)
for suffix in ['', '.step0.csv','.step1.csv','.budget.csv']:
 name='result.csv'+suffix;assert (baseline/name).read_bytes()==(enabled/name).read_bytes(),name
assert not list(baseline.glob('*.mol_*.csv'))
for name,extra,error,ranks in [('count',' mol_trace_count=17,\n','Invalid MOL trace cell count',0),('node',' mol_trace_count=1,mol_trace_nodes(1)=3,\n','Invalid MOL trace node',0),('cell',' mol_trace_count=1,mol_trace_cells(:,1)=5,1,1,\n','Invalid MOL trace cell',0),('mpi',trace,'MOL cell trace currently requires serial gas',2)]:run(name,extra,ranks,error)
def rows(path):
 with path.open() as f:out=list(csv.DictReader(f))
 assert out and all(None not in row and all(math.isfinite(float(v)) for v in row.values()) for row in out)
 return out
stage=rows(enabled/'result.csv.mol_stage.csv');sync=rows(enabled/'result.csv.mol_sync.csv')
for row in stage:
 if float(row['aperture'])>0:
  assert all(float(row[k])>0 for k in ['rho_left','rho_right','sound_left','sound_right']),row
stencil=rows(enabled/'result.csv.mol_stencil.csv')
stencil_lookup={}
for row in stencil:
 key=tuple(row[k] for k in ['event','probe','axis','side','offset'])
 assert key not in stencil_lookup
 stencil_lookup[key]=row
 if float(row['rho'])>0:
  assert float(row['p'])>0 and float(row['sound'])>0,row
  assert abs(math.fsum(float(v) for k,v in row.items() if k.startswith('Y_'))-1)<2e-12,row
for row in stage:
 key=tuple(row[k] for k in ['event','probe','axis','side'])
 assert all(key+(str(o),) in stencil_lookup for o in [-2,-1,0,1])
 if float(row['aperture'])>0:
  for offset,name in [(-1,'sound_left'),(0,'sound_right')]:
   assert float(stencil_lookup[key+(str(offset),)]['sound'])==float(row[name])
 # The selected cell is the right donor at its lower face and left at upper.
 center=stencil_lookup[key+(str(-int(row['side'])),)]
 assert float(center['p'])==float(row['pressure_cell_Pa'])
groups=collections.defaultdict(list)
for r in stage:groups[(r['event'],r['probe'],r['component'])].append(r)
max_residual=0
for key,rs in groups.items():
 assert len(rs)==6 and {(int(r['axis']),int(r['side'])) for r in rs}=={(a,s) for a in [1,2,3] for s in [0,1]}
 r=rs[0];dx=.25/(2 if int(r['node'])==2 else 1)
 face_rate=math.fsum((1 if int(x['side'])==0 else -1)*float(x['aperture'])*float(x['centroid_hydro_flux']) for x in rs)/(dx*float(r['vfrac']))
 expected=face_rate+float(r['wall_rate']);actual=float(r['hydro_rate']);scale=max(1,abs(expected),abs(actual))
 residual=abs(actual-expected)/scale;max_residual=max(max_residual,residual);assert residual<2e-12,(key,residual)
assert {int(r['component']) for r in stage}==set(range(1,6+len([x for x in rows(enabled/'result.csv')[0] if x.startswith('rhoY_')])))
assert {int(r['phase']) for r in sync}=={0,1}
reacting_off=run('reacting_disabled',reacting=True)
reacting_on=run('reacting_enabled',trace,reacting=True)
for suffix in ['', '.step0.csv','.step1.csv','.budget.csv']:
 name='result.csv'+suffix
 assert (reacting_off/name).read_bytes()==(reacting_on/name).read_bytes(),('reacting',name)
reaction=rows(reacting_on/'result.csv.mol_reaction.csv');max_reaction_residual=0
for row in reaction:
 old=float(row['old_state']);new=float(row['updated_state']);dt=float(row['interval_s'])
 expected=dt*(float(row['forcing_rate'])+float(row['reaction_rate']))
 residual=abs(new-old-expected)/max(1,abs(old),abs(new),abs(expected))
 max_reaction_residual=max(max_reaction_residual,residual)
 assert residual<2e-12,row
assert {int(r['stage']) for r in reaction}=={1,2}
(w/'checks.json').write_text(json.dumps({'status':'passed','qualification':False,'runs':runs,'reaction_rows':len(reaction),'maximum_scaled_reaction_balance_residual':max_reaction_residual,'stencil_rows':len(stencil),'stage_rows':len(stage),'sync_rows':len(sync),'balance_groups':len(groups),'maximum_scaled_hydro_balance_residual':max_residual},indent=2)+'\n')
print(w)

"""Exercise explicit EB gas state admission and owner-local distribution."""
from pathlib import Path
import argparse
import subprocess
import os
import csv
import json
import hashlib
import tempfile

parser = argparse.ArgumentParser()
for option in ('serial', 'mpi', 'mpiexec', 'work'):
    parser.add_argument('--' + option, required=True)
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
Path(args.work).mkdir(parents=True, exist_ok=True)
w = Path(tempfile.mkdtemp(prefix='initial_state_', dir=Path(args.work).resolve()))
serial = Path(args.serial).resolve()
mpi = Path(args.mpi).resolve()
env = os.environ.copy()
base="""&eb_tree
 root_shape=2,2,2,child_count=0,
 dual_energy_enabled=.true.,dual_energy_eta2=2,allow_negative_energy=.false.,
 time_integrator='sdc',mol_mode=.true.,sdc_iterations=2,
 hydro_method='plm_ctu',riemann_solver='pelec',boundary_faces=6*'outflow',
 regular_density=.31,cut_density=.35,initial_pressure=250000,
 geometry_type='plane',plane_axis='z',plane_position=.31,
 time_step=1.e-9,final_time=1.e-9,maximum_steps=1,
 chemistry_enabled=.false.,transport_enabled=.false.,
 composition_count=3,composition_species='H2','O2','N2',composition_mole_fractions=2,1,3.76,
 output_interval_steps=1,
"""
results=[]
def run(name,extra='',ranks=0,expected=None,template=None):
 d=w/name;d.mkdir(exist_ok=False)
 inp=d/'input.nml';inp.write_text((base if template is None else template)+extra+f" output_file='{d}/result.csv',\n/\n")
 cmd=([str(serial)] if not ranks else [args.mpiexec,'-n',str(ranks),str(mpi)])+[str(inp)]
 with (d/'run.log').open('w') as log:p=subprocess.run(cmd,env=env,cwd=d,stdout=log,stderr=subprocess.STDOUT,timeout=60)
 results.append({'name':name,'command':cmd,'returncode':p.returncode})
 (w/'progress.json').write_text(json.dumps(results,indent=2)+'\n')
 if expected:
  assert p.returncode!=0 and expected in (d/'run.log').read_text(),name
  assert not (d/'result.csv').exists(),name
 else:assert p.returncode==0,(name,(d/'run.log').read_text())
 return d
seed=run('seed')
rows=list(csv.DictReader((seed/'result.csv.step0.csv').open()))
names=[s[5:] for s in rows[0] if s.startswith('rhoY_')]
rows.sort(key=lambda r:tuple(int(r[k]) for k in ['k','j','i']))
lines=['PELEF_EB_TREE_STATE_V1',f'1 {len(names)}',*names,'1 2 2 2']
for r in rows:
 vals=[r[x] for x in ['rho','rhou','rhov','rhow','rhoE']]+[r['rhoY_'+n] for n in names]+[r['rho_internal_energy_J_m3'],r['temperature']]
 lines.append(' '.join([r[x] for x in ['node','i','j','k']]+vals))
assert len(rows)==8
state=w/'state.dat';state.write_text('\n'.join(lines)+'\n')
for name,ranks in [('serial',0),('mpi2',2),('mpi3',3)]:
 d=run(name,f" initial_state_file='{state}',\n",ranks)
 out=list(csv.DictReader((d/'result.csv').open()));ref=list(csv.DictReader((seed/'result.csv').open()))
 assert len(out)==len(ref)
 for a,b in zip(out,ref):
  for k in a:
   assert abs(float(a[k])-float(b[k]))<=1.e-12*max(1.,abs(float(b[k]))),(name,k,a[k],b[k])

def change_cell(lines, column, value):
 fields=lines[-1].split();fields[column]=value;lines[-1]=' '.join(fields)

for name,mutate,msg in [
 ('species',lambda x:x.__setitem__(2,'WRONG'),'Initial state species order mismatch'),
 ('missing',lambda x:x.pop(),'Missing initial cell'),
 ('order',lambda x:x.__setitem__(-1,x[-2]),'Initial cell ordering mismatch'),
 ('mass',lambda x:change_cell(x,4,'0.7'),'Initial species mass mismatch'),
 ('temperature',lambda x:change_cell(x,-1,'1000'),'Initial temperature mismatch'),
 ('nonfinite',lambda x:change_cell(x,4,'NaN'),'Nonfinite initial cell'),
 ('extra',lambda x:x.append('extra'),'Unexpected initial state trailing record')]:
 ls=lines.copy();mutate(ls);f=w/(name+'.dat');f.write_text('\n'.join(ls)+'\n')
 for ranks in [0,2]:run(name+str(ranks),f" initial_state_file='{f}',\n",ranks,msg)
# A child on a distinct owner imports covered coarse donors and all fine
# cells; compare both the initial snapshot and one evolved step across ranks.
multi_lines=lines.copy();multi_lines[1]=f'2 {len(names)}'
multi_lines.append('2 2 2 2')
seed_values=lines[-8].split()[4:]
for k in range(1,3):
 for j in range(1,3):
  for i in range(1,3):
   multi_lines.append(' '.join(['2',str(i),str(j),str(k)]+seed_values))
multi=w/'multi.dat';multi.write_text('\n'.join(multi_lines)+'\n')
template=base.replace('child_count=0,', 'child_count=1,parent=1,ratio=2,lower(:,1)=1,1,1,upper(:,1)=1,1,1,')
reference=None
for ranks in (0,2,3):
 d=run('multi'+str(ranks),f" initial_state_file='{multi}',\n",ranks,template=template)
 if reference is None:
  reference=d
 else:
  for file in ('result.csv.step0.csv','result.csv'):
   assert (reference/file).read_bytes()==(d/file).read_bytes(),(ranks,file)
# Restriction of opposite fine velocities leaves kinetic variance in the
# covered parent, even though its mean momentum is zero.
variance_lines=multi_lines.copy()
for index,line in enumerate(variance_lines):
 fields=line.split()
 if len(fields)!=4+5+len(names)+2:continue
 node,i,j,k=map(int,fields[:4])
 if node==2 or (node,i,j,k)==(1,1,1,1):
  rho=float(fields[4]);fields[8]=str(float(fields[8])+0.5*rho*4.0)
  if node==2:fields[5]=str(rho*(2.0 if i==1 else -2.0))
  variance_lines[index]=' '.join(fields)
variance=w/'variance.dat';variance.write_text('\n'.join(variance_lines)+'\n')
for ranks in (0,2,3):
 run('variance'+str(ranks),f" initial_state_file='{variance}',\n",ranks,template=template)
# Covered donor temperature is an independently restricted EOS seed.
restricted=variance_lines.copy()
for index,line in enumerate(restricted):
 fields=line.split()
 if len(fields)==4+5+len(names)+2 and tuple(map(int,fields[:4]))==(1,1,1,1):
  fields[-1]=str(float(fields[-1])+1.0);restricted[index]=' '.join(fields)
path=w/'restricted_temperature.dat';path.write_text('\n'.join(restricted)+'\n')
for ranks in (0,2):
 run('restricted_temperature'+str(ranks),f" initial_state_file='{path}',\n",ranks,template=template)
for label,node,indices,delta,message in [
 ('covered_deficit',1,(1,1,1),-1.0,'Initial covered energy deficit'),
 ('leaf_excess',1,(2,2,2),1.0,'Initial held energy mismatch')]:
 modified=multi_lines.copy()
 for index,line in enumerate(modified):
  fields=line.split()
  if len(fields)==4+5+len(names)+2 and tuple(map(int,fields[:4]))==(node,*indices):
   fields[8]=str(float(fields[8])+delta);modified[index]=' '.join(fields)
 path=w/(label+'.dat');path.write_text('\n'.join(modified)+'\n')
 for ranks in (0,2):
  run(label+str(ranks),f" initial_state_file='{path}',\n",ranks,message,template=template)
for ranks in (0,2):
 run('restart_conflict'+str(ranks),f" initial_state_file='{state}',restart_file='unused.chk',\n",ranks,
     'Initial state file conflicts with restart')
(w/'checks.json').write_text(json.dumps({'qualification':False,'status':'passed','results':results,'binary_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (serial,mpi)},'source_sha256':hashlib.sha256((root/'app/pelef_reactive_eb_tree_3d.F90').read_bytes()).hexdigest()},indent=2)+'\n')
print('PASS initial-state serial/MPI admission smoke',len(results))

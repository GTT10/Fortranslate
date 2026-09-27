"""Compare pinned PeleC pc_initdata with the SI gamma-law Fortran application.

Only the initialization kernel is compiled, with minimal AMReX storage adapters.
This is not a PeleC time-integration comparison.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

PELEC = 'bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8'
PHYSICS = '6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b'


def function(text, name):
    start = text.index('void\n' + name + '(')
    brace = text.index('{', start)
    depth = 1
    end = brace + 1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}')
        end += 1
    return text[start:end]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('reference', 'exe', 'work'):
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--cxx', default='/usr/bin/g++')
    args = parser.parse_args()
    args.exe = args.exe.resolve(strict=True)
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='hit_reference_', dir=args.work.resolve()))
    records = {}

    def blob(repo, commit, path, sha):
        data = subprocess.check_output(['git', '-C', str(repo), 'show', commit + ':' + path])
        actual = hashlib.sha256(data).hexdigest()
        assert actual == sha, (path, actual)
        (work / Path(path).name).write_bytes(data)
        records[path] = dict(commit=commit, sha256=actual)
        return data.decode()

    prob = blob(args.reference, PELEC, 'Exec/RegTests/HIT/prob.H',
                '54400bea0fedb634bdd0cfc61653b4276d71f6e85015c7bfdf1dd90ad373a4a4')
    utility = blob(args.reference / 'Submodules/PelePhysics', PHYSICS,
                   'Source/Utility/Utilities/Utilities.H',
                   '7d19e0003be3e94627db9778f42cb9ba857cebd5a3f555b678c74a9f00733ebe')
    blob(args.reference, PELEC, 'Exec/RegTests/HIT/hit_ic_4_32.dat',
         '983ebe8d4db570e56a634346871993e346ba6b3bb7eebf02ca8f39a6e5194776')
    harness = r'''
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
namespace amrex {
using Real = double;
struct GeometryData {
 double lower[3]={0,0,0}, spacing[3];
 const Real* ProbLo() const { return lower; }
 const Real* CellSize() const { return spacing; }
};
template<class T> struct Array4 {
 T* data;
 T& operator()(int i,int j,int k,int n) const { return data[(((k*4)+j)*4+i)*8+n]; }
};
}
constexpr int NUM_SPECIES=1;
enum {URHO,UEINT,UTEMP,UFS,UMX,UMY,UMZ,UEDEN};
struct ProbParmDevice {
 double rho0,eint0,T0,Linput,forcing_u0=0,forcing_v0=0,forcing_w0=0;
 int inres=32;
 double *d_xarray,*d_xdiff,*d_uinput,*d_vinput,*d_winput;
};
namespace pele::physics::utilities {
''' + function(utility, 'locate') + '\n}\n' + function(prob, 'pc_initdata') + r'''
int main() {
 std::ifstream input("hit_ic_4_32.dat");
 std::string line;
 std::getline(input,line);
 std::vector<double> knots(32),width(32),u(32768),v(32768),w(32768);
 const double R=8.31446261815324e7/28.97, gamma=1.4, T=300;
 const double rms=.1*std::sqrt(gamma*R*T)/std::sqrt(3.);
 for(int n=0;n<32768;n++) {
   if(!std::getline(input,line)) return 2;
   std::replace(line.begin(),line.end(),',',' ');
   std::istringstream row(line);
   double x,y,z;
   if(!(row>>x>>y>>z>>u[n]>>v[n]>>w[n])) return 3;
   if(n<32) knots[n]=x;
   u[n]*=rms; v[n]*=rms; w[n]*=rms;
 }
 for(int n=1;n<32;n++) width[n]=knots[n]-knots[n-1];
 width[0]=width[1];
 ProbParmDevice p;
 p.rho0=1.013e6/(R*T); p.eint0=R*T/(gamma-1); p.T0=T;
 p.Linput=knots[31]+.5*width[31];
 p.d_xarray=knots.data(); p.d_xdiff=width.data();
 p.d_uinput=u.data(); p.d_vinput=v.data(); p.d_winput=w.data();
 amrex::GeometryData geometry;
 for(int n=0;n<3;n++) geometry.spacing[n]=6.283185307179586/4;
 std::vector<double> values(4*4*4*8);
 amrex::Array4<double> state{values.data()};
 std::cout<<"rho,rhou,rhov,rhow,rhoE,temperature\n"<<std::setprecision(17);
 for(int k=0;k<4;k++) for(int j=0;j<4;j++) for(int i=0;i<4;i++) {
   pc_initdata(i,j,k,state,geometry,p);
   std::cout<<state(i,j,k,URHO)*1000<<','<<state(i,j,k,UMX)*10<<','
    <<state(i,j,k,UMY)*10<<','<<state(i,j,k,UMZ)*10<<','
    <<state(i,j,k,UEDEN)*.1<<','<<state(i,j,k,UTEMP)<<'\n';
 }
}
'''
    (work / 'reference.cpp').write_text(harness)
    env = dict(os.environ)
    for key in ('LD_LIBRARY_PATH', 'LIBRARY_PATH', 'DISPLAY'):
        env.pop(key, None)

    def run(command, log):
        result = subprocess.run(command, cwd=work, env=env, capture_output=True,
                                text=True, timeout=120)
        (work / log).write_text(result.stdout + result.stderr)
        assert result.returncode == 0, result.stdout + result.stderr
        return result.stdout

    run([args.cxx, '-std=c++17', '-O0', '-o', 'reference', 'reference.cpp'], 'compile.log')
    run([str(work / 'reference')], 'reference.csv')
    (work / 'input.nml').write_text('''&reactive_3d
nx=4,ny=4,nz=4
x_upper=.06283185307179586,y_upper=.06283185307179586,z_upper=.06283185307179586
final_time=1e-20,maximum_steps=2,time_integrator='sdc',sdc_iterations=2
problem='uniform_reactor',thermo_model='gamma_law'
reconstruction='characteristic_ppm',riemann_solver='pelec',boundary_condition='periodic'
chemistry_enabled=.false.,transport_enabled=.false.
species_diffusion_enabled=.false.,barodiffusion_enabled=.false.
initial_temperature=300,initial_pressure=101300
initial_velocity_x=0,initial_velocity_y=0,initial_velocity_z=0
output_file='fortran.csv'
/
&gamma_law gamma=1.4 /
&reactive_hit input_file='hit_ic_4_32.dat',inres=32,coordinate_scale=.01,mach_t0=.1,uin_norm=1 /
''')
    run([str(args.exe.resolve()), str(work / 'input.nml')], 'fortran.log')
    with (work / 'reference.csv').open() as stream:
        reference = list(csv.DictReader(stream))
    with (work / 'fortran.csv').open() as stream:
        actual = list(csv.DictReader(stream))
    assert len(reference) == len(actual) == 64
    for index, row in enumerate(actual):
        # The C++ adapter publishes x fastest, then y, then z.
        for key, cell in zip(('x', 'y', 'z'), (index % 4, index // 4 % 4, index // 16)):
            expected = (cell + .5) * .06283185307179586 / 4
            assert abs(float(row[key]) - expected) < 1e-16, row
        assert float(row['time']) == 1e-20
    rho = 101300 / (8314.46261815324 / 28.97 * 300)
    momentum = rho * .1 * math.sqrt(1.4 * 8314.46261815324 / 28.97 * 300) / math.sqrt(3)
    scales = dict(rho=rho, rhou=momentum, rhov=momentum, rhow=momentum,
                  rhoE=101300/.4, temperature=300)
    errors = {key: max(abs(float(a[key])-float(b[key]))/scale
                       for a, b in zip(actual, reference)) for key, scale in scales.items()}
    assert all(math.isfinite(e) and e < 1e-12 for e in errors.values()), errors
    records['normalized_errors'] = errors
    records['fortran_executable_sha256'] = hashlib.sha256(args.exe.read_bytes()).hexdigest()
    records['compiler'] = subprocess.check_output([args.cxx, '--version'], text=True).splitlines()[0]
    records['comparison'] = 'pc_initdata at t=0 versus Fortran after 1e-20 s, transport and chemistry disabled'
    (work / 'summary.json').write_text(json.dumps(records, indent=2) + '\n')
    print('PASS: pinned HIT initialization, 64 cells, SI normalized errors:', errors)
    print(work)


if __name__ == '__main__':
    main()

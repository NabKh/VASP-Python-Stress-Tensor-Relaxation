#!/usr/bin/python

#SBATCH -N 4
#SBATCH -J name
#SBATCH --ntasks-per-node=32
#SBATCH -t 24:00:00
#SBATCH -p genoa


'''
This script is a utility to iteratively perform structure relaxation in VASP until a desired stress tensor is achieved.

Instructions:
1. Define the pressure (stress) tensor in the script (Setpress variable).
2. Define the estimated elastic modulus (E) and Poisson's ratio (v) in the script.
3. Set ISIF=2 in the INCAR file.
4. Copy POSCAR to poscar.0 for initial conditions.
5. Configure mpiexec command according to your VASP directory and job management system.
6. Change directory to your job directory and submit your job using 'qsub fixpressure.py'.

Suggestions for accurate force calculation:
- Use a high ECUT and k-points density. Refer to the VASP manual for details.
- Use PREC=High or Accurate for precision settings.

Suggestions for computational efficiency:
- Use LWAVE=.TRUE. and the default ISTART and ICHARG settings to reuse wavefunctions.

General algorithm:
1. Run VASP with ISIF=2 to obtain the current pressure tensor.
2. Use generalized Hooke's law with the input elastic modulus to estimate the required POSCAR modifications.
3. Repeat until the pressure difference meets the convergence criteria.

Contact for bugs or issues: n.khossossi@tudelft.nl
'''

import sys,os
import subprocess
import shutil
import string
import numpy as np
import linecache

#######################You need to configurate the following variables##############
Setpress= np.array([120,0,0,0,0,0])
#Set the pressure tensor (xx yy zz xy yz zx) in kB. Note pressure=-stress in VASP
presscirt= 0.1                        
#convergency criteria for pressure, unit in kB
E=2790                           
#Young's modulus in kB
v=0.21                          
#Possion ratio
G=E/(2+2*v)
#Shear modulus in kB
imax=30
#maximum iteration cycles
mpiexec='mpirun -n `cat $PBS_NODEFILE | wc -l` vasp_std > results.txt'
#directory of your VASP program 
os.chdir(os.getenv('PBS_O_WORKDIR', ''))
#these PBS command may depend on your job management system.
##################################################################################


########do not alter the following codes unless you know what you are doing#######
shutil.copy ('poscar.0','POSCAR')
subprocess.call("echo 'starting calculation' > pressure.all", shell=True)
subprocess.call(mpiexec, shell=True)
#initial run to calculate pressure tensor
subprocess.call("cat OSZICAR > oszicar.all", shell=True)
subprocess.call("cat OUTCAR > outcar.all", shell=True)
subprocess.call("grep 'Total CPU time used (sec):' OUTCAR | tail -n 1 > end.txt", shell=True)

iteration=0
while iteration<=imax:
 iteration=iteration+1
 P_coeff=1-iteration*1.0/imax
 #graduately reduce P_coeff to avoid flucations in convergency
 print  'iteration=',iteration
 subprocess.call("grep 'Total CPU time used (sec):' OUTCAR | tail -n 1 > end.txt", shell=True)
 notempty=os.path.getsize('end.txt')
 os.remove('end.txt')
 #scan if the previous calculation is completed
 if notempty:
  subprocess.call(" grep 'in kB' OUTCAR | tail -n 1 > pressure.txt", shell=True)
  subprocess.call(" grep 'in kB' OUTCAR | tail -n 1 >> pressure.all", shell=True)
  for line in open("pressure.txt"):  
   press=np.array([ float(x) for x in line.split()[2:]])
   print  'pressure=',press
  os.remove('pressure.txt')
  #read pressure from OUTCAR
  addpress=Setpress-press[0:6]
  print  'adding pressure:',addpress
  #calculate the additional pressure needed to achieve Setpress
  abs_P=max([abs(y) for y in addpress])
  if (abs_P<presscirt):
   subprocess.call("echo 'aborting calculation for the convergence is reached' >> pressure.all", shell=True)
   break
   #stop if the additional pressure is too small
  else:
   #copy CONTCAR to POSCAR, 
   shutil.copy ('CONTCAR','POSCAR')
   POS=linecache.getlines('POSCAR')
   M=np.zeros((3,3))
   addM=np.zeros((3,3))
   M[0,:]=np.array([float(x) for x in POS[2].split()])
   M[1,:]=np.array([float(x) for x in POS[3].split()])
   M[2,:]=np.array([float(x) for x in POS[4].split()])
   #read current POSCAR matrix
   addM[0,0]=(addpress[0]-v*(addpress[1]+addpress[2]))/(-E)
   addM[1,1]=(addpress[1]-v*(addpress[0]+addpress[2]))/(-E)
   addM[2,2]=(addpress[2]-v*(addpress[1]+addpress[0]))/(-E)
   addM[0,1]=addM[1,0]=addpress[3]/(-G)/2
   addM[1,2]=addM[2,1]=addpress[4]/(-G)/2
   addM[0,2]=addM[2,0]=addpress[5]/(-G)/2
   #calculate how much strain we need to add according to addpress
   #beware that VASP defines compressive pressure as positive, hence /E becomes /(-E) and /G becomes /(-G)
   addM=np.diag([1,1,1])+addM*P_coeff
   M=np.dot(M,addM)
   print  'adjusting matrix to:'
   print  M
   np.savetxt("M.txt",M)
   linecache.clearcache()
   addPOS=linecache.getlines('M.txt')
   POS[2:5]=addPOS
   subprocess.call("rm M.txt", shell=True)
   fo = open("POSCAR", "w+")
   fo.writelines(POS)
   fo.close()
   posname='poscar.'+str(iteration)
   shutil.copy ('POSCAR',posname)

   #run vasp to calculate pressure for this POSCAR
   subprocess.call(mpiexec, shell=True)
   subprocess.call("cat OSZICAR >> oszicar.all", shell=True)
   #run VASP program, configure the shell command according to your linux environment
 else:
  print 'calculation aborted'
  #previous vasp simulation is not completed, calculation abort
subprocess.call("rm WAVECAR", shell=True)
subprocess.call("rm outcar.all", shell=True)

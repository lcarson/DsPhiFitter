import os,time,sys, pwd, random, datetime, shutil
todaysdate = datetime.date.today()

def makeSubmitScript(toyDir, seed, n, wall , brstring):
    initDir=os.getcwd();
    dirname=toyDir
    os.chdir(dirname)
    os.system("echo '#!/bin/bash' >> BatchSubmit_"+seed+".sh")
    os.system("echo '#PBS -l cput="+wall+":00:00' >> BatchSubmit_"+seed+".sh")
    os.system("echo '#PBS -l walltime="+wall+":00:00' >> BatchSubmit_"+seed+".sh")
    os.system("echo 'source /data/lhcb/sw/scripts/lbsetup-cvmfs.sh' >> BatchSubmit_"+seed+".sh")
    os.system("echo '. SetupProject.sh Gaudi ROOT' >> BatchSubmit_"+seed+".sh")
    os.system("echo 'cd "+dirname+"/../' >> BatchSubmit_"+seed+".sh")
    if n=='1':
        os.system("echo 'bin/run -H -m Ds2PhiPi:Ds2KKPi:Ds2PiPiPi:Ds2KPiPi -d s21:s21r1:s24:s26 -S -M -c DsD0:DsPhi:DsPhiSide -X -Z "+brstring+" -s "+seed+"'>> BatchSubmit_"+seed+".sh")
    else:  
        os.system("echo 'bin/run -H -m Ds2PhiPi:Ds2KKPi:Ds2PiPiPi:Ds2KPiPi -d s21:s21r1:s24:s26 -S -M -c DsD0:DsPhi:DsPhiSide -Z "+brstring+"-t "+n+" -T toysDir/ -s "+seed+"'>> BatchSubmit_"+seed+".sh")

    os.chdir(initDir)

def makePlotScript(toyDir, wall, brstring):
    initDir=os.getcwd();
    dirname=toyDir
    os.chdir(dirname)
    os.system("echo '#!/bin/bash' >> MakePlots.sh");
    os.system("echo '#PBS -l cput="+wall+":00:00' >> MakePlots.sh");
    os.system("echo '#PBS -l walltime="+wall+":00:00' >> MakePlots.sh");
    os.system("echo 'source /data/lhcb/sw/scripts/lbsetup-cvmfs.sh' >> MakePlots.sh");
    os.system("echo '. SetupProject.sh Gaudi ROOT' >> MakePlots.sh");

    os.system("echo 'cd "+dirname+"/../' >> MakePlots.sh")
    os.system("echo 'bin/run -H -m Ds2PhiPi:Ds2KKPi:Ds2PiPiPi:Ds2KPiPi -d s21:s21r1:s24:s26 -S -M -c DsD0:DsPhi:DsPhiSide -Z "+brstring+"-t -T toysDir/ -X'>> MakePlots.sh");
    os.chdir(initDir)

basedir = "/data/lhcb/users/hadavizadeh/B2DsPhi/ToyStudies"

desc = raw_input("Enter a description of the toys you are running (today's date will automatically be included): ")

if desc == '':
    print 'You entered a blank string so I am using "toys"'
    desc = 'toys'

toyDir = basedir + '/' + str(todaysdate) + "_" + desc
if os.path.exists(toyDir):
    print
    print '*** WARNING: directory', toyDir, 'already exists! ***'
    decision = raw_input('*** Input Y to overwrite, anything else to abort: ').lower()
    if decision != 'y':
        sys.exit()
    else:
        shutil.rmtree(toyDir)
        print 'Directory has been removed and will be recreated'
print 'Creating', toyDir
os.mkdir(toyDir)
print
print
brtype = str(raw_input("Fit yields, Single branching fraction, or four branching fractions? type y, s, or f: "))
if brtype != 'y' and brtype != 'f' and brtype != 's':
    print 'Decision not reccognised: '+ brtype
brstring = ''
if brtype=='s':
    brstring = ' -R '
if brtype=='f':
    brstring = ' -r '

mergedInput = str(raw_input("Merged plots? y or other: "))
if mergedInput=='y':
    brstring += ' -D '

# Get seed
seedInput = int(raw_input("Give first random seed: "))
# Get Number of batch jobs
numBatchJobs = int(raw_input("Num batch jobs: "))
# Get Number of Toys per job
numMCStudyToysPerBatchJob = int(raw_input("Num MCStudyToys per batch job (>1): "))
# Get wall time for job
numWalltime = int(raw_input("How many hours should it run for? "))

print 'Copying the code into the new toy directory.' 
# Make a copy of the code there ready to run - then don't have to wait for each toy to complete
os.system('mkdir '+toyDir+'/Bu2DsPhi_Fitter_copy')
os.system('mkdir '+toyDir+'/Bu2DsPhi_Fitter_copy/src')
os.system('mkdir '+toyDir+'/Bu2DsPhi_Fitter_copy/scripts')
os.system('mkdir '+toyDir+'/Bu2DsPhi_Fitter_copy/bin')
os.system('mkdir '+toyDir+'/Bu2DsPhi_Fitter_copy/pdfs')
os.system('mkdir '+toyDir+'/Bu2DsPhi_Fitter_copy/results')
os.system('mkdir '+toyDir+'/Bu2DsPhi_Fitter_copy/toysDir')
os.system('mkdir '+toyDir+'/Bu2DsPhi_Fitter_copy/toysDir/plots')
os.system('mkdir '+toyDir+'/Bu2DsPhi_Fitter_copy/toysDir/gen_vals')
os.system('mkdir '+toyDir+'/Bu2DsPhi_Fitter_copy/batch')
os.system("cp -r /home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/src/* "+toyDir+"/Bu2DsPhi_Fitter_copy/src")
os.system("cp -r /home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/scripts/* "+toyDir+"/Bu2DsPhi_Fitter_copy/scripts")
os.system("cp -r /home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/bin/* "+toyDir+"/Bu2DsPhi_Fitter_copy/bin")
os.system("cp -r /home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/pdfs/* "+toyDir+"/Bu2DsPhi_Fitter_copy/pdfs")


jobCounter=1
toyDirName = toyDir+"/Bu2DsPhi_Fitter_copy"
os.system('cd '+toyDirName + '/batch')

for toySeed in range(seedInput,numBatchJobs*numMCStudyToysPerBatchJob+seedInput,numMCStudyToysPerBatchJob):
    makeSubmitScript(toyDirName+ '/batch', str(toySeed), str(numMCStudyToysPerBatchJob),str(numWalltime), brstring) 
    
    print 'Submitting job '+str(jobCounter)+" with seed "+str(toySeed)
    time.sleep(0.1)
    os.chdir(toyDirName + '/batch')
    os.system('qsub -N Toy_'+str(toySeed)+' BatchSubmit_'+str(toySeed)+'.sh')
    print
    jobCounter+=1

makePlotScript(toyDirName+'/batch',str(numWalltime), brstring)
print 
print "Toys can be found in:"
os.system('pwd')

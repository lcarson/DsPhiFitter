import os,time,sys, pwd, random, datetime, shutil
todaysdate = datetime.date.today()

def makeSubmitScript(sysDir, seed, N, wall, brstring,inputMode):
    initDir=os.getcwd();
    dirname=sysDir
    os.chdir(dirname)
    os.system("echo '#!/bin/bash' >> BatchSubmit_"+seed+".sh");
    os.system("echo '#PBS -l cput="+wall+":00:00' >> BatchSubmit_"+seed+".sh");
    os.system("echo '#PBS -l walltime="+wall+":00:00' >> BatchSubmit_"+seed+".sh");
    os.system("echo 'source /data/lhcb/sw/scripts/lbsetup-cvmfs.sh' >> BatchSubmit_"+seed+".sh");
    os.system("echo '. SetupProject.sh Gaudi ROOT' >> BatchSubmit_"+seed+".sh");
    os.system("echo 'cd "+dirname+"/../' >> BatchSubmit_"+seed+".sh")
    if inputMode=='DsPhi':
        os.system("echo 'bin/run -m Ds2PhiPi:Ds2KKPi:Ds2PiPiPi:Ds2KPiPi -S -M -c DsD0:DsPhi:DsPhiSide -Z -Q -H -X -f "+N+":"+brstring+" -s "+seed+"'>> BatchSubmit_"+seed+".sh");
    elif inputMode=='DKst0':
        os.system("echo 'bin/run -m D2PiKPi -S -M -c DD0:DKst0:DKst0Side -Z -Q -H -X -f "+N+":"+brstring+" -s "+seed+"'>> BatchSubmit_"+seed+".sh");
    os.chdir(initDir)


basedir = "/data/lhcb/users/hadavizadeh/B2DsPhi/Sensitivity"

desc = raw_input("Enter a description of the fits you are running (today's date will automatically be included): ")

if desc == '':
    print 'You entered a blank string so I am using "sense"'
    desc = 'sense'

sysDir = basedir + '/' + str(todaysdate) + "_" + desc
if os.path.exists(sysDir):
    print
    print '*** WARNING: directory', sysDir, 'already exists! ***'
    decision = raw_input('*** Input Y to overwrite, anything else to abort: ').lower()
    if decision != 'y':
        sys.exit()
    else:
        shutil.rmtree(sysDir)
        print 'Directory has been removed and will be recreated'
print 'Creating', sysDir
os.mkdir(sysDir)
print


# Get seed
inputMode     = str(raw_input("DsPhi or DKst0? "))
if inputMode!='DsPhi' and inputMode!= 'DKst0':
    print "Mode not recognised: " + inputMode
    print "Assuming you meant DsPhi"

# Get seed
seedInput     = int(raw_input("Give first random seed: "))
# Get Number of batch jobs
numToysPerJob = int(raw_input("Num toys per job: "))
# Get wall time for job
numJobsPerBr  = int(raw_input("Num jobs per BR value: "))
# Get Number toys per job


#brValues = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10]
#brValues = [0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0]
brValues = [0.0,3.0,6.0,9.0,12.0,15.0,18.0,21.0]
#brValues = [0.1,0.2,0.4,0.6,0.8]
numBrValues = len(brValues)
print "Number of BR values used: " + str(numBrValues)
print brValues
print

numBatchJobs  = numBrValues*numJobsPerBr
print "Total number of batch jobs: " + str(numBatchJobs)

numWalltime   = int(raw_input("How many hours should each job run for? "))

print 'Copying the code into the new toy directory.' 
# Make a copy of the code there ready to run - then don't have to wait for each toy to complete
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/src')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/scripts')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/bin')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/pdfs')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/sensitivityDir')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/sensitivityDir/gen_vals')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/sensitivityDir/text_vals')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/batch')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/roodatasets')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/data')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/results')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/data/FullSel')
os.system("cp -r /home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/src/* "+sysDir+"/Bu2DsPhi_Fitter_copy/src")
os.system("cp -r /home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/scripts/* "+sysDir+"/Bu2DsPhi_Fitter_copy/scripts")
os.system("cp -r /home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/bin/* "+sysDir+"/Bu2DsPhi_Fitter_copy/bin")
os.system("cp -r /home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/pdfs/* "+sysDir+"/Bu2DsPhi_Fitter_copy/pdfs")


jobCounter=1
sysDirName = sysDir+"/Bu2DsPhi_Fitter_copy"
os.system('cd '+sysDirName + '/batch')

toySeed = seedInput
os.chdir(sysDirName + '/batch')

for br in brValues:
    print "===> Branching ratio value: " + str(br)

    for i in range(0,numJobsPerBr):
        makeSubmitScript(sysDirName+'/batch', str(toySeed), str(numToysPerJob), str(numWalltime), str(br),inputMode)
        print 'Submitting job '+str(jobCounter)+" with seed "+str(toySeed)
        time.sleep(0.1)
        os.system('qsub -N Sensitivity_'+str(toySeed)+' BatchSubmit_'+str(toySeed)+'.sh')
        #os.system('echo -N Sensitivity_'+str(toySeed)+' BatchSubmit_'+str(toySeed)+'.sh')
        print
        jobCounter+=1
        toySeed+=1

print 
print "Toys can be found in:"
os.system('pwd')

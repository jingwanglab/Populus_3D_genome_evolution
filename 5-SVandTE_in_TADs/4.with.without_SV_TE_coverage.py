import pandas as pd

## 1. prepare data file
spes=['Plas','Psim','Ppse','Pwua','Psze','Pyun','Prot','Pdav','Pqio','Pade']
for spe in spes:
    d1 = pd.read_csv('./'+spe+'/'+spe+'.con.boundary.sv.coverage.csv',sep='\t',header=None)
    d1['sv'] = d1.apply(lambda x: 'sv-' if x[4]==0 else 'sv+', axis=1)
    d1['sv'].value_counts()
    d1 = d1[[0,1,2,'sv']]
    
    d2 = pd.read_csv('./'+spe+'/'+spe+'.high.con.boundary.sv.coverage.csv',sep='\t',header=None)
    d2['sv'] = d2.apply(lambda x: 'sv-' if x[4]==0 else 'sv+', axis=1)
    d2['sv'].value_counts()
    d2 = d2[[0,1,2,'sv']]
    
    d3 = pd.read_csv('./'+spe+'/'+spe+'.spe.boundary.sv.coverage.csv',sep='\t',header=None)
    d3['sv'] = d3.apply(lambda x: 'sv-' if x[4]==0 else 'sv+', axis=1)
    d3['sv'].value_counts()
    d3 = d3[[0,1,2,'sv']]
    
    d1.to_csv('../TE_SV/'+spe+'/'+spe+'.con.boundary.sv.bed',sep='\t',header=None,index=False)
    d2.to_csv('../TE_SV/'+spe+'/'+spe+'.high.con.boundary.sv.bed',sep='\t',header=None,index=False)
    d3.to_csv('../TE_SV/'+spe+'/'+spe+'.spe.boundary.sv.bed',sep='\t',header=None,index=False)



tes=['all_TE.bed','TE_DNA.bed','TE_helitron.bed','TE_LTR.bed']
cmd=open('TE_class.sv.cmd','w')
for spe in spes:
    for te in tes:
        t = te.split('.')[0]
        cmd.write('bedtools coverage -nonamecheck -b /usr_storage2/stt/tad/TE/'+spe+'/'+te+' -a /usr_storage2/stt/tad/TE_SV/'+spe+'/'+spe+'.spe.boundary.sv.bed > /usr_storage2/stt/tad/TE_SV/'+spe+'/'+'spe.boundary.sv.'+t+'.coverage'+'\n')
        cmd.write('bedtools coverage -nonamecheck -b /usr_storage2/stt/tad/TE/'+spe+'/'+te+' -a /usr_storage2/stt/tad/TE_SV/'+spe+'/'+spe+'.con.boundary.sv.bed > /usr_storage2/stt/tad/TE_SV/'+spe+'/'+'con.boundary.sv.'+t+'.coverage'+'\n')
        cmd.write('bedtools coverage -nonamecheck -b /usr_storage2/stt/tad/TE/'+spe+'/'+te+' -a /usr_storage2/stt/tad/TE_SV/'+spe+'/'+spe+'.high.con.boundary.sv.bed > /usr_storage2/stt/tad/TE_SV/'+spe+'/'+'high.con.boundary.sv.'+t+'.coverage'+'\n')

## 2. run cmd file


## 3. combine results        
spes=['Plas','Psim','Ppse','Pwua','Psze','Pyun','Prot','Pdav','Pqio','Pade']
tes=['all_TE','TE_DNA','TE_helitron','TE_LTR']
for te in tes:
    df = pd.DataFrame(columns=[0,1,2,3,4,5,6,7,'species'])
    for spe in spe:
        d1= pd.read_csv('/usr_storage2/stt/tad/TE_SV/'+spe+'/'+'high.con.boundary.sv.'+te+'.coverage',sep='\t',header=None)
        d2= pd.read_csv('/usr_storage2/stt/tad/TE_SV/'+spe+'/'+'con.boundary.sv.'+te+'.coverage',sep='\t',header=None)
        d3= pd.read_csv('/usr_storage2/stt/tad/TE_SV/'+spe+'/'+'spe.boundary.sv.'+te+'.coverage',sep='\t',header=None)
        d1['type']='H-conserved'
        d2['type']='Conserved'
        d3['type']='Diverged'
        d = pd.concat([d1,d2,d3])
        d['species']=spe
        df = pd.concat([df,d])
        
    df.to_csv('all_spe_bound_SV.'+te+'.coverage.csv',index=False)
        
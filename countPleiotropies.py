# Author: Juan A. Rodriguez
# Creation Date:
# Last changed at Time-stamp: <2019-02-06 16:53:27 (jrodriguez)>
## Description:

'''
Checks the risk haplotype in the haplotype freq files
coming from plink output in script haplotypes_pleiotropy.sh
'''
import sys

def readFreqHapFile(path):
    '''
    '''
    hapFreqs={}
    with open(path, 'r') as fhtp:
        for line in fhtp:
            (key, val) = line.split(' ')
            hapFreqs[(key)] = float(val.strip())
    return hapFreqs

def checkRiskHaploFreqs(hapFreqs, riskHap, pleio, agon, antagon):
    '''
    '''
    if len(hapFreqs) == 2:
#        print '2 Haplotypes'
        if riskHap in hapFreqs:
#            print 'AGON'
            agon.append(pleio)
        else:
#            print 'ANT'
            antagon.append(pleio)
    else:
#        print len(hapFreqs), 'Haplotypes'
        justFreqs = sorted(hapFreqs.values())[:-2]
        minim = sum(justFreqs)
        if minim > 0.1:
            #print(hapFreqs)
            sys.exit(['This is exceeding the frequency threshold.'])
        elif riskHap in hapFreqs and hapFreqs[riskHap] > 0.1:
#            print 'AGON'
            agon.append(pleio)
        else:
#            print 'ANT'
            antagon.append(pleio)

#He eliminado el filtro de edad debido a que las infecciosas
#tienen una edad inferior a la esperada (1 año por defecto)
def writeFiles(agon, antagon, path_out):
    '''
    '''
    print(path_out)
    if agon:
        with open(path_out + "Agon", 'a') as agF:
            for i in agon:
                agF.write(str(i))
                agF.write('\n')
    if antagon:
        with open(path_out + "Antagon", 'a') as antF:
            for i in antagon:
                antF.write(str(i))
                antF.write('\n')

def main():
    """
    main function
    """
    agon=[]
    antagon=[]
#     path = '/home/jrodriguez/Projects/PLEIOTROPY/Hapls_Output/rs4766578_rs3184504.fhtp'
    path = sys.argv[1]
#    riskHap='AC' ## ${riskHap}
    riskHap = sys.argv[2]
#     car='rs1333042	Coronary heart disease	A	rs6475606	Intracranial aneurysm	T	LATE	EARLY'
    pleio_rel = list(sys.argv[3:-1]) #pleiotropic relation from haplotypes_pleiotropies.sh ${CAR}

    out_path = sys.argv[-1]
## Input of pleiotropic relation from haplotypes_pleiotropies.sh ${CAR}
    #age = sys.argv[4] #age threshold
# LD Threshold/file analyzed name
    #ld = sys.argv[5]
    hapFreqs = readFreqHapFile(path)
    checkRiskHaploFreqs(hapFreqs, riskHap, pleio_rel, agon, antagon)
    writeFiles(agon, antagon, out_path)

if __name__ == "__main__":
    exit(main())

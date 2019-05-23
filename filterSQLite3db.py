#!/usr/bin/python3
# Author: Juan A. Rodriguez
# Creation Date: jue nov 7 16:09:38 CET 2013
# Last changed at Time-stamp: <2014-10-09 12:59:46 (jrodriguez)>
## Description:

'''
Refilter the pleiotropic table uploaded to the MySQL
(allPleiotropies)
'''

import sqlite3
import collections
import sys
import itertools

DBpath = '../results/GWASpleiotropies.sqlite'
r2=0.8

"""
Create Iterative object
"""

def connectPleiotropy(path, chr, r2):
    conn = sqlite3.connect(path) # Make aconnection object9cursorObject = conn.cursor () # Create a cursor  object
    pleio = conn.cursor() # Create a cursor  object
    statement = 'SELECT * FROM filteredPairs WHERE R2>={0} AND CHR={1}'.format(r2,chr)
    pleio.execute(statement)
    return pleio

'''
Makes the pairs of the diseases in the pleiotropy.
'''

def makePairsDictionary(pleio):
    dis={}
    for row in pleio:
        parejas = (sorted([row[2],row[8]]))
        name_disPair = '_'.join(parejas)
        if name_disPair in dis.keys():
            dis[name_disPair].append(row[0])
        else:
            dis[name_disPair] = []
            dis[name_disPair].append(row[0])
    return dis


"""
Function to make a dictionary with the repeated
events. Con esta funcion tendremos las señales duplicadas
en un diccionario; con el id y la posicion. Necesitamos saber
cual es el r2 mas elevado y nos
quedaremos con ese ID.
Makes a default dictionary to find the repeats.
"""

def makeRepeatsDictionary(dis, pleio, path, chr, r2, quitar):
    reps={}
    for k, v in dis.items():
        #print(k, v)
        if len(v) > 1:
            reps.update({k:v})
    for k, v in reps.items():
        #print('Repetidos',k ,v)
        if len(v) > 2:
            sql = "SELECT POS1, ID, R2 FROM filteredPairs WHERE ID IN (%s)"
            in_p = ', '.join(list(map(lambda x: '%s', v)))
            sql = sql % (in_p)
            v = tuple(v)
            pleio.execute(sql % v)
            RepPleio = list(pleio.fetchall())
            multiSigns = [(pair[0], pair[1]) for pair in itertools.combinations(RepPleio, 2)]
            #print(multiSigns)
            for i in multiSigns:
                if abs(i[0][0]-i[1][0]) < 200000:
                    #print('Same signal', i)
                    ir2 = tuple(sorted(i, key=lambda item: item[2]))[:-1][0][1]#coge el que tiene el R2 mas bajo para eliminar
                    #print(ir2)
                    quitar.add(ir2)
                    #print('\n Se eliminan:', ir2, '\n')
        else:
            sits=[]
            dit={}
            pleio = connectPleiotropy(path, chr, r2)
            for p in pleio:
                if p[0] in v:
                    pos=min(p[6],p[12])
                    clau=p[0]
                    ld=p[13]
                    l=[pos,ld]
                    dit.update({clau:l})
            #print(dit)
            for t,y in dit.items():
                sits.append(y)
            #print(sits)
            dif=abs(sits[0][0]-sits[1][0])
            if dif < 200000:
                #print(dit.keys(), 'misma senhal')
                eliminamos=sorted(dit.items(), key=lambda e: e[1][1])[:-1]#coge el que tiene el R2 mas bajo para eliminar
                quitar.add(eliminamos[0][0])
                #print('\n Se eliminan:', eliminamos, '\n')
    return quitar

#def keep_only_one_pleiotropy_diseases():
    #Iterate through each pair of diseases
    #pleio.execute(statement) statement select * from table where disease1 = tal and disease2 = tal otro
    #for row in subtable of pair of diseases:
        # for next row :
            # get info from positions
            # if min - max otro o min otro menos max 200,000 nt
                # mirar el R2 mas grande de los dos
                # si el que tiene el R2 mas pequeño esta en la lista
                    #eliminar de la lista la fila y guardar el mayor
                # si no es de esta manera
                    #guardar la fila del que tiene el R2 Mas grande


    #AL FINAL DEL PROCESO TENDRAS EL SUBSET DE ESA PAREJA DE ENFERMEDADES CON LAS filas
    # DE LOS SNPS CON R2 MAS ALTO Y  ESTAN A distancia considerable ENTRE ELLOS

"""
MAIN Script
"""

def main():
    """
    main function
    """
    DBpath = '../results/GWASpleiotropies.sqlite'
    r2=0.8
    quitar=set()
    for chr in range(1,23):
        pleiotropies = connectPleiotropy(DBpath, chr, r2)
        diseases = makePairsDictionary(pleiotropies)
        quitar = makeRepeatsDictionary(diseases, pleiotropies, DBpath, chr, r2, quitar)
        #print quitar
    q = tuple(quitar)
    if len(q) != 0:
        conn = sqlite3.connect(DBpath)
        curs = conn.cursor()
        for element in q:
            sql = '''DELETE FROM filteredPairs WHERE ID = {0}'''.format(element)
            curs.execute(sql)
        conn.commit()
        conn.close()
    else:
        print('No hay repetidos')


if __name__ == "__main__":
    exit(main())

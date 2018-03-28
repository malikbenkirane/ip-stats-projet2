# -*- coding: utf-8 -*-

def read_fasta(filename):
    """
    Lis le fichier fileneame (au format fasta)
    Retourne la liste des sequences lues

    Arguments:
    filename -- chemin du fichier a lire (string)
    """
    with open(filename, 'r') as fd:
        seqs = []
        for line in fd:
            if line[0] == '>':
                continue
            seqs.append(line[:-1])
    return seqs

def read_distances(filename):
    """
    Lis le fichier de distance : 3 colonnes
    * ID de la première position
    * ID de la seconde position
    * distance
    Retourne liste des triplets lus

    Arguments:
    filename -- chemin du fichier a lire (string)
    """
    with open(filename, 'r') as fd:
        distances = {}
        for line in fd:
            t = tuple(line.split())
            id_1 = int(t[0])
            id_2 = int(t[1])
            dist = float(t[2])
            distances[(id_1, id_2)] = dist
    return distances

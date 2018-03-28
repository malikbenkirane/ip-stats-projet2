# -*- coding: utf-8 -*-

import numpy as np

A = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
     'T', 'V', 'W', 'Y', '-']


def estimate_cooc(data):
    """
    Calcul les co-occurences puis les poids associés.
    Retourne le dictionnaire des poids indexé par les couples possibles d'acides
    aminés, a chaque index correspond un dictionnaire des poids associées aux
    co-occurences indexées par les positions.

    Arguments:
    data -- liste de sequences (cf. module data, read_fasta)
    """

    assert len(np.unique([len(seq) for seq in data])) == 1

    # initialisation
    L = len(data[0])
    co_occurences = {
        (a, b): {} for a in A for b in A
    }

    #comptage
    for sequence in data:
        for i, a in enumerate(sequence):
            for j, b in enumerate(sequence):
                if not (i, j) in co_occurences[(a, b)]:
                    co_occurences[(a, b)][(i, j)] = 0
                co_occurences[(a, b)][(i, j)] += 1

    # normalisation
    q, M = len(A), len(data)
    for coaa in co_occurences:
        for copos in co_occurences[coaa]:
            co_occurences[coaa][copos] = \
                (co_occurences[coaa][copos] + 1 / q) / (M+q)

    return co_occurences

def get_cooc(a, b, i, j, M, cooc):
    """
    Retourne le poids associé aux co-occurences pour a, b aux positions i, j
    Simple adjustement parce que cooc est représenté par un dictionnaire.
    Il s'agit de traiter les couples qui n'apparaissent durant le comptage.

    Arguments:
    a, b -- acide aminés (comme décrits par A)
    i, j -- positions (entiers)
    M  -- nombre de séquences (pour normaliser)
    cooc -- matrice co-occurences (comme dans estimate_cooc)
    """

    q = len(A)
    if (i, j) not in cooc[(a,b)]:
        return 1 / (q * (M+q))
    return cooc[(a, b)][(i, j)]


def mutual_information(cooc, pswm, M):
    """
    Calcul et retourne la matrice d'information mutuelle.
    Nécessite d'avoir précalculé les poids associés aux co-occurences
    et la matrice de poids spécifique (PSWM).

    Arguments:
    cooc -- poids des co-occurences (voir estimate_cooc)
    pswm -- matrice PSWM (voir pswm.estimate_pswm)
    M -- nombre de séquences
    """

    def e(a, b, i, j):
        c = get_cooc(a, b, i, j, M, cooc)
        return c * (np.log2(c) - np.log2(pswm[a][i]) - np.log2(pswm[b][j]))

    L = len(pswm[A[0]])
    mim = np.zeros((L,L))

    for i in range(L):
        for j in range(L):
            mim[i,j] = sum(e(a, b, i, j) for a, b in cooc)

    return mim


def induced_contacts(mim, distances):
    """
    Determine les fractions des paires en contact en fonction du nombre de paires
    sélectionnées par quantité d'information mutuelle décroissante. Deux paires
    sont considérées en contact si leur distance est inférieur a 8.
    Retourne la liste des fractions du nombre de paire en contact.

    Arguments:
    mim -- matrice d'information mutuelle (voir mutual_information)
    distances -- voir data.read_distances
    """

    L = len(mim)
    imsort = [(i//L, i%L) for i in reversed(np.argsort(mim, axis=None))]
    
    def warn(i, j):
        warnings += 1
        print('Warning (%i): position ignorée, donnée manquante (%i, %i)' \
              % (warnings, i, j))

    ncontacts, fracs = 0, []
    warnings = 0
    for n, (i, j) in enumerate(imsort):
        d = 0
        # lecture de la distance
        if (i, j) in distances:
            d = distances[(i, j)]
        else: 
            if i != j:
                if (j, i) in distances:
                    d = distances[(j, i)]
                else:
                    warn(i, j)
                    continue
        # test contact
        if d < 8:
            ncontacts += 1
        # ajout de la fraction
        fracs.append(ncontacts / (n + 1))

    print('-- Warnings : %i positions ignorées (%.2f%%)' \
          %(warnings, warnings*100/n))

    return fracs

def real4():
    """
    Implémentation de la réalisation "quatrième fonction" partie co-évolution
    sur les données de data/Dtrain.txt avec les distances data/distances.txt
    """
    from data import read_fasta, read_distances
    from pswm import estimate_pswm
    from matplotlib import pyplot as plt

    # lecture des données
    Dtrain = read_fasta('data/Dtrain.txt')
    distances = read_distances('data/distances.txt')

    # préparation
    P = estimate_pswm(Dtrain)
    CO = estimate_cooc(Dtrain)
    MIM = mutual_information(CO, P, len(Dtrain))
    # calcul des fractions
    fracs = induced_contacts(MIM, distances)

    # plot
    fig, ax = plt.subplots()
    ax.plot(fracs)
    ax.set_xlabel('quantité de corrélation')
    ax.set_ylabel('fraction contacts')
    fig.show()
    return fig

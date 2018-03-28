# -*- coding: utf-8 -*-

import numpy as np

A = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
     'T', 'V', 'W', 'Y', '-']

def estimate_pswm(data):
    """
    Calcul et retourne la matrice de poids specifique de position
    La matrice est represente par un dictionnaire dont les cles sont
    les elements de A (acide aminées). Le poids de chaque acide aminé
    a la colonne i est accessible dans la liste contenue a sa clé.

    (pswm stands for position-specific weight matrix)

    Arguments:
    data -- matrice de sequences d'acides amines
    """

    assert len(np.unique([len(seq) for seq in data])) == 1

    # compteur d'occurences initial
    cols = len(data[0])
    occurences = {}
    for aa in A:
        occurences[aa] = [0] * cols

    # print(occurences.keys())

    # comptage
    for sequence in data:
        for col, aa in enumerate(sequence):
            occurences[aa][col] += 1

    # calcul pwm
    pswm = {}
    M = len(data)
    q = len(A)
    for aa in A:
        pswm[aa] = [(na + 1) / (M + q) for na in occurences[aa]]

    return pswm


def shannon_pswm(i, pswm):
    """
    Calcul l'entropie relative d'une position i dans la séquence et l'acide aminé
    correspondant.
    Retourne l'acide aminé pour cette position et l'entropie correspondante
    sous la forme d'unc couple

    Arguments:
    i -- position évaluée
    pswm -- matrice PSWM telle qu'estimée par estiate_pswm
    """

    q = len(A)
    shannon = np.log2(q) + sum([pswm[aa][i] * np.log2(pswm[aa][i]) for aa in A])
    aamax = A[np.argsort([pswm[aa][i] for aa in A])[-1]]

    return aamax, shannon


def real2():
    """
    Implémentation de la réalisation "deuxième fonction" partie PSWM
    sur les données de data/Dtrain.txt
    """
    from data import read_fasta
    from matplotlib import pyplot as plt

    Dtrain = read_fasta('data/Dtrain.txt')
    P = estimate_pswm(Dtrain)
    S = [shannon_pswm(i, P) for i in range(len(Dtrain[0]))]
    ai, Si = list(zip(*S))
    top3 = np.argsort(Si)[-3:]

    fig, ax = plt.subplots()
    ax.set_xlabel('Position')
    ax.set_ylabel('Entropie Relative')
    ax.set_xticks(range(0,48,2))
    ax.plot(Si)
    for i in top3:
        ax.plot((i, i), (5, Si[i]), linestyle='dashed', color='grey')
        ax.plot((i, i), (0, 1), linestyle='dashed', color='grey')
    ax.set_ylim(bottom=0)
    fig.show()
    return fig


def null_model(pswm):
    """
    Retourne le dictionnaire dont les clés sont les acides aminé de A
    et dont les valeurs sont les composantes du modèle nulle.

    Arguments:
    pswm -- matrice PSWM telle qu'estimé par estimate_pswm
    """
    L = len(pswm[A[0]])
    return {
        aa: sum(pswm[aa]) / L
        for aa in A
    }


def log_odds(sequence, pswm, f0):
    """
    Calcul et retourne la log-vraissemblance d'une séquence a partir des 
    composantes du modèle nulle et de la matrice  pswm de la famille a comparer

    Arguments:
    sequence -- chaine de caractères qui représente la séquence
    pswm -- matrice de poids spécifique calculable avec estimate_pswm
    f0 -- dictionnaire des composantes du modèle nulle, elles sont calculable 
          avec la fonction null_model
    """
    return np.sum(
        np.log2(pswm[aa][i]) - np.log2(f0[aa])
        for i, aa in enumerate(sequence)
    )


def sliding_odds(sequences, pswm):
    """
    Calcul les log-vraissemblance de sequences longue en faisant glisser
    une fenêtre glissante sur ces séquences.
    Retourne la liste des log-vraissemblances calculés.

    Arguments:
    sequences -- liste de séquence a comparer avec Dtrain
    pswm -- matrice de poids spécifiques telle que décrite par estimate_pswm
    """
    f0 = null_model(pswm)
    L, windows, odds = len(pswm[A[0]]), [], []
    for sequence in sequences:
        # data/testseq.txt ne contient qu'une seule sequence
        # on ne passera qu'une seule fois dans cette boucle
        N = len(sequence)
        windows += [sequence[shift:shift+L] for shift in range(N)]
    for window in windows:
        odds.append(log_odds(window, pswm , f0))
    return odds


def real4():
    """
    Implémentation de la réalisation "quatrième fonction" partie PSWM
    avec testseq.txt et Dtrain.txt
    """
    from data import read_fasta
    from matplotlib import pyplot as plt
    P = estimate_pswm(read_fasta('data/Dtrain.txt'))
    odds = sliding_odds(read_fasta('data/test_seq.txt'), P)
    fig, ax = plt.subplots()
    ax.plot(odds)
    ax.set_ylabel('log-odds')
    ax.set_xlabel('window shift')
    ax.plot([0,len(odds)], [0,0], color='grey', linestyle='dashed')
    ax.set_xticks(range(0,len(odds),10))
    x,y = odds.index(max(odds)), max(odds)
    plt.tight_layout()
    fig.show()
    return fig

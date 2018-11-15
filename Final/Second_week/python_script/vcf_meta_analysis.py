#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import argparse
import os

RED = "#EB6231"
BLUE = "#5D80B4"
GREEN = "#8FB03E"


def build_array_from_vcf(path, vcf_name, _cutoff):
    vcf_file = open("{0}/{1}".format(path, vcf_name), 'r')
    variant_list = []
    nbr_indiv = 0
    for vcf_line in vcf_file:
        # Pour chaque SNP, faire la boucle: regarde SNP par SNP
        if vcf_line[0] == "#":
            if vcf_line[1] != "#":
                nbr_indiv = len(vcf_line.split("\t")) - 9
                print("{0} individuals found in the vcf file.".format(nbr_indiv))
        else:
            line_list = vcf_line.replace('\n', '').split("\t")[9:9 + nbr_indiv]

            sparse_line = np.zeros(nbr_indiv, dtype=np.int8)
            for i, genotype in enumerate(line_list):
                # Pour chaque individu, regarder son génotype au SNP étudié
                value = genotype.count("1")
                # Compte le nb de 1: 0 (0/0 homozygote référence),1 (0/1 ou 1/0 hétérozygote), 2 (1/1 homozygote varaint)
                if value > 0:
                    # Ajoute le génotype à value et si value sup à 0, conserver la valeur de l'individu et l'ajouter à la somme de tous les individus
                    sparse_line[i] = value
            if args.c >= np.sum(sparse_line) > 0:
                # Si la somme pour le SNP est comprise entre le cut-off et 0, garder le SNP
                variant_list.append(sparse_line)
    vcf_file.close()

    return np.array(variant_list, dtype=np.int8)


def epistasis(array):
    va = np.sum([np.var(col) for col in array])
    s2 = np.var([np.sum(row) for row in array.T])
    return s2 / va


def bootstrap(nbr_snps, array, bootstrap_resample=1000):
    # nbr_snps: nombre snp stop
    # array: comme définit précédemment dans vcf_analysis
    # bootstrap: nombre de tirages de nbr_snps

    ligne, colonne = array.shape
    liste_s2_Va = []
    for i in range(bootstrap_resample):
        choice = np.random.choice(range(ligne), nbr_snps, replace=False)
        # Tirage aléatoire dans tous les snps présents dans array (syn ou nonsyn), pour le nombre de snps stop
        sub_array = array[choice, :]
        # Sous-tableau aléatoire
        assert ((sub_array[0, :] == array[choice[0], :]).all())
        # Vérifier que la première ligne du subarray correspond bien au premier tirage aléatoire
        liste_s2_Va.append(epistasis(sub_array))

    # return: une liste de s2/Va
    return np.array(liste_s2_Va)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--stop', required=True, type=str,
                        dest="stop", metavar="<vcf>",
                        help="The relative name of the stop .vcf file")
    parser.add_argument('-s', '--syn', required=True, type=str,
                        dest="syn", metavar="<vcf>",
                        help="The relative name of the synonymous .vcf file")
    parser.add_argument('-n', '--nonsyn', required=True, type=str,
                        dest="non_syn", metavar="<vcf>",
                        help="The relative name of the non-synonymous .vcf file")
    parser.add_argument('-c', '--cutoff', required=False, default=1, type=int,
                        dest="c", metavar="<panel>",
                        help="The cut-off for minor allele count")

    args = parser.parse_args()

    stop_array = build_array_from_vcf(os.getcwd(), args.stop, args.c)
    print("Computed stop variants array for the whole population.")
    syn_array = build_array_from_vcf(os.getcwd(), args.syn, args.c)
    print("Computed synonymous variants array for the whole population.")
    non_syn_array = build_array_from_vcf(os.getcwd(), args.non_syn,  args.c)
    print("Computed non-synonymous variants array for the whole population.")

    tsv_file = open("{0}/meta_analysis_{1}.tsv".format(os.getcwd(), args.c), 'w')
    tsv_file.write('\t'.join(['Population', 'NbrIndividuals', 'NbrStops', 's^2/Va', 'Pvalue', 'Significant']) + '\n')

    nbr_indiv = stop_array.shape[1]
    print("{0} individuals".format(nbr_indiv))

    stop_epistasis = epistasis(stop_array)

    nbr_stops = stop_array.shape[0]
    if nbr_stops > 0:
        syn_bootstrap = bootstrap(nbr_stops, syn_array)
        non_syn_bootstrap = bootstrap(nbr_stops, non_syn_array)

        my_dpi = 128
        fig = plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
        p_val = len([1 for epi in syn_bootstrap if epi < stop_epistasis]) / syn_bootstrap.shape[0]

        result_line = [str(nbr_indiv), str(nbr_stops), str(stop_epistasis), str(p_val)]
        if p_val < 0.05:
            result_line += ["True"]
            print("\t p-value={0:3g}, the statistical test is significant ".format(p_val))
        else:
            result_line += ["False"]
            print("\t p-value={0:3g}, the statistical test is not significant ".format(p_val))
        tsv_file.write('\t'.join(result_line) + '\n')

        plt.title('$p={0:3g}$'.format(p_val))
        bins = 25
        syn_hist, _, _ = plt.hist(syn_bootstrap, bins, density=1, facecolor=BLUE, alpha=0.4, label='Syn')
        non_syn_hist, _, _ = plt.hist(non_syn_bootstrap, bins, density=1, facecolor=GREEN, alpha=0.4, label='NonSyn')
        y_max = 1.2 * max((max(syn_hist), max(non_syn_hist)))
        plt.ylim((0, y_max))
        plt.plot((stop_epistasis, stop_epistasis), (0, y_max), label="LoF", linewidth=3, color=RED)
        plt.xlabel('$\sigma^2/V_{A}$')
        plt.ylabel('density')
        plt.legend()
        plt.tight_layout()
        plt.savefig("{0}/{1}.analysis.{2}.svg".format(os.getcwd(), args.stop, args.c), format="svg")
        plt.savefig("{0}/{1}.analysis.{2}.png".format(os.getcwd(), args.stop, args.c), format="png")
        plt.close()
    else:
        print('\tNo stop variants for this population')

    tsv_file.close()
    print("Analysis completed")

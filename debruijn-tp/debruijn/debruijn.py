#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import matplotlib
import networkx as nx
import operator
import random
random.seed(9001)
import statistics

__author__ = "Romain Pholoppe"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Romain Pholoppe"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Romain Pholoppe"
__email__ = "pholoppero@eisti.eu"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    fichier = open(fastq_file, "rt")
    for lignes in fichier:
        yield next(fichier).strip("\n")
        next(fichier)
        next(fichier)
    fichier.close()


def cut_kmer(read, kmer_size):
    for i in range(len(read) + 1 - kmer_size):
        kmer = read[i:i + kmer_size]
        yield(kmer)


def build_kmer_dict(fastq_file, kmer_size):
    dico=dict()
    for ligne in read_fastq(fastq_file):
        temp = cut_kmer(ligne, kmer_size)
        for kmer in temp:
            if kmer in dico:
                dico[kmer] += 1
            else :
                dico[kmer] = 1
    return dico


def build_graph(kmer_dict):
    graph = nx.DiGraph()
    for kmer in kmer_dict:
        precedant = kmer[:-1]
        suivant = kmer[1:]
        graph.add_edge(precedant, suivant, weight = kmer_dict[kmer])
    return graph



def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass

def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    noeuds = []
    for n in graph.nodes() :
        if not graph.pred[n]:
            noeuds.append(n)
    return noeuds

def get_sink_nodes(graph):
    noeuds = []
    for n in graph.nodes() :
        if not graph.succ[n]:
            noeuds.append(n)
    return noeuds

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for S_noeud in starting_nodes:
        for E_noeud  in ending_nodes:
            for path in nx.all_simple_paths(graph,S_noeud,E_noeud):
                temp = []
                for n in range(len(path)) :
                    temp.append(path[n][0])
                temp.append(path[n][1])
                contigs.append(("".join(temp), len(temp)))
    return contigs


def fill(text, width=80):
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contigs_list, output_file):
    save = open(output_file,"w+")
    for i in range(len(contigs_list)):
        save.write(">contig_%d len=%d\n"%(i,contigs_list[i][1]))
        save.write(fill(contigs_list[i][0],80))
        save.write("\n")
    save.close()





#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

if __name__ == '__main__':
    main()

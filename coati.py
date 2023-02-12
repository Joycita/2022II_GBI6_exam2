from Bio import Phylo
from Bio import SeqIO
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Entrez
from Bio import SeqIO
from Bio import GenBank 
import csv 
import re 
import matplotlib
import matplotlib.pyplot as plt
from Bio.Align.Applications import ClustalwCommandline
import os

##Primera función.- 
def fasta_downloader():
    '''
    La siguiente función con el nombre fasta_downloader nos permite ejecutar:
    Primero con el with open el archivo coati.txt que se encuentra dentro de la carpeta data las vamos abrir como coatlist y posteior a esto vamos a convertirla en lista, luego se hace uso del entrez de Biopythony con la
    finalidad de obtener la secuencia de cada identificador en el formato genbank y finalmente guardamos los datos en un documento .gb en la carpeta data
    '''
    with open('data\coati.txt', 'r+') as coatlist:
        id_coati= coatlist.readlines()
        coati= []
        Entrez.email= "jozzfer14@gmail.com"
        for b in id_coati:
            handle = Entrez.efetch(db="nucleotide", id=b, rettype="gb", retmode="text") 
            seq_record= SeqIO.read(handle, "genbank")
            coati.append(seq_record)
            SeqIO.write(coati, "data\coati.gb", "genbank")
            
##Segunda función
from Bio.Align.Applications import ClustalwCommandline
import os
def alignment():
    '''
    La siguiente función con el nombre alignment nos permite ejecutar:
     Primero transformar el archivo que se tiene de tipo genbank en formato fasta definiendo las variables que se presentan a continuación y posterior a esto se realizó un alineamiento de las secuencias empleando ClustalW con ayuda de la variable culstalw_exe y así poder acceder al programa ClustalW2, una vez realizado el alineamiento el resultado es guardado en los archivos coati.aln y coati.dnd.
     '''
    sequences= SeqIO.parse('data\coati.gb', 'genbank')
    SeqIO.write(sequences, "data\coati.fasta", "fasta")  
    cline = ClustalwCommandline("clustalw2", infile="data/coati.fasta")
    align = AlignIO.read("data/coati.aln", "clustal")
    print(align)
    
##Tercera función
def tree():
    """
    La función con el nombre tree es capaz de realizar un cálculo de distancias haciendo uso del archivo aln para poder formar un árbol filogenético y su resultado poder guardarlo en un archivo pdf, para esto se toman encuenta las siguiente varibles como alignment, distance_calculation, dmatrix...
    """
    with open("data\coati.aln", "r") as f:
        alignment  = AlignIO.read(f,"clustal")
        distance_calculation= DistanceCalculator("identity")
        dmatrix= distance_calculation.get_distance(alignment)
        tree_constructor=DistanceTreeConstructor(distance_calculation)
        tree=tree_constructor.build_tree(alignment)
        tree.rooted= True
        Phylo.draw(tree)
        plt.show()
        plt.savefig("data/coati_phylotree.pdf")

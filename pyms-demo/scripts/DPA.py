#!/usr/bin/env python
# coding: utf-8

# ## Example: Within-state alignment of peak lists from multiple experiments
#
# In this example the experiments ``a0806_077``, ``a0806_078``, and ``a0806_079``
# prepared in the previous example will be aligned, and therefore the notebook
# ``Multiple_Experiments.ipynb`` must be run first to create the files
# ``a0806_077.expr``, ``a0806_078.expr``, ``a0806_079.expr``. These files contain
# the post-processed peak lists from the three experiments.
#
# First, determine the directory to the experiment files and import the required functions.

# In[18]:


import pathlib
output_directory = pathlib.Path(".").resolve() / "output"

from pyms.DPA.PairwiseAlignment import PairwiseAlignment, align_with_tree
from pyms.DPA.Alignment import exprl2alignment
from pyms.Experiment import load_expr


# Define the input experiments list.

# In[19]:


exprA_codes = ["a0806_077", "a0806_078", "a0806_079"]


# Read the experiment files from disk and create a list of the loaded |Experiment| objects.

# In[20]:


expr_list = []

for expr_code in exprA_codes:
    file_name = output_directory / "experiments" / f"{expr_code}.expr"
    expr = load_expr(file_name)
    expr_list.append(expr)


# Define the within-state alignment parameters.

# In[21]:


Dw = 2.5  # rt modulation [s]
Gw = 0.30 # gap penalty


# Convert each |Experiment| object is converted into an |Alignment| object with
# the function |exprl2alignment()|.

# In[22]:


F1 = exprl2alignment(expr_list)


# In this example, there is only one experimental condition so the alignment
# object is only for within group alignment (this special case is called
# 1-alignment). The variable ``F1`` is a Python list containing three alignment
# objects.
#
# Perform pairwise alignment. The class |pyms.DPA.Class.PairwiseAlignment|
# calculates the similarity between all peaks in one sample with those of another sample.
# This is done for all possible pairwise alignments (2-alignments).

# In[23]:


T1 = PairwiseAlignment(F1, Dw, Gw)


# The parameters for the alignment by dynamic programming are: ``Dw``, the
# retention time modulation in seconds; and ``Gw``, the gap penalty. These
# parameters are explained in detail in [1]_.
#
# The output of |PairwiseAlignment| (``T1``) is an object which contains the
# dendrogram tree that maps the similarity relationship between the input
# 1-alignments, and also 1-alignments themselves.
#
# The function |align_with_tree()| then takes the object ``T1`` and aligns the
# individual alignment objects according to the guide tree.

# In[24]:


A1 = align_with_tree(T1, min_peaks=2)


# In this example, the individual alignments are three 1-alignments, and the
# function |align_with_tree()| first creates a 2-alignment from the two most
# similar 1-alignments and then adds the third 1-alignment to this to create
# a 3-alignment.
#
# The parameter ``min_peaks=2`` specifies that any peak column of the data
# matrix that has fewer than two peaks in the final alignment will be dropped.
# This is useful to clean up the data matrix of accidental peaks that are not
# truly observed over the set of replicates.
#
# Finally, the resulting 3-alignment is saved by writing alignment tables
# containing peak retention times (``rt.csv``) and the corresponding peak areas
# (``area.csv``). These are plain ASCII files in CSV format.

# In[25]:


A1.write_csv(
		output_directory / "within_state_alignment" / 'a_rt.csv',
		output_directory / "within_state_alignment" / 'a_area.csv',
		)


# The file ``area1.csv`` contains the data matrix where the corresponding peaks are aligned in the columns and each row corresponds to an experiment.
# The file ``rt1.csv`` is useful for manually inspecting the alignment.
#
# ## Example: Between-state alignment of peak lists from multiple experiments
#
# In the previous example the list of peaks were aligned within a single
# experiment with multiple replicates ("within-state alignment"). In practice, it
# is of more interest to compare the two experimental states.
#
# In a typical experimental setup there can be multiple replicate experiments on
# each experimental state or condition. To analyze the results of such an
# experiment statistically, the list of peaks need to be aligned within each
# experimental state and also between the states. The result of such an alignment
# would be the data matrix of integrated peak areas. The data matrix contains a
# row for each sample and the number of columns is determined by the number of
# unique peaks (metabolites) detected in all the experiments.
#
# In principle, all experiments could be aligned across conditions and replicates
# in the one process. However, a more robust approach is to first align
# experiments within each set of replicates (within-state alignment), and then to
# align the resulting alignments (between-state alignment) [1]_.
#
# This example demonstrates how the peak lists from two cell states are aligned.
#
# * Cell state A, consisting of three aligned experiments
# (``a0806_077``, ``a0806_078``, and ``a0806_079``), and
# * Cell state B, consisting of three aligned experiments
# (``a0806_140``, ``a0806_141``, and ``a0806_142``).
#
# These experiments were created in the notebook ``Multiple_Experiments.ipynb``.
#
# First, perform within-state alignment for cell state B.

# In[26]:


exprB_codes = ["a0806_140", "a0806_141", "a0806_142"]

expr_list = []

for expr_code in exprB_codes:
    file_name = output_directory / "experiments" / f"{expr_code}.expr"
    expr = load_expr(file_name)
    expr_list.append(expr)

F2 = exprl2alignment(expr_list)
T2 = PairwiseAlignment(F2, Dw, Gw)
A2 = align_with_tree(T2, min_peaks=2)

A2.write_csv(
		output_directory / "within_state_alignment" / 'b_rt.csv',
		output_directory / "within_state_alignment" / 'b_area.csv',
		)


# ``A1`` and ``A2`` are the results of the within group alignments for cell state A and B, respectively.
# The between-state alignment can be performed as follows alignment commands:

# In[27]:


# Define the within-state alignment parameters.
Db = 10.0 # rt modulation
Gb = 0.30 # gap penalty

T9 = PairwiseAlignment([A1,A2], Db, Gb)
A9 = align_with_tree(T9)

A9.write_csv(
		output_directory / "between_state_alignment" / 'rt.csv',
		output_directory / "between_state_alignment" / 'area.csv')


# Store the aligned peaks to disk.

# In[28]:


from pyms.Peak.List.IO import store_peaks

aligned_peaks = A9.aligned_peaks()
store_peaks(aligned_peaks, output_directory / "between_state_alignment" / 'peaks.bin')


# In this example the retention time tolerance for between-state alignment is
# greater compared to the retention time tolerance for the within-state alignment
# as we expect less fidelity in retention times between them. The same functions
# are used for the within-state and between-state alignment. The result of the
# alignment is saved to a file as the area and retention time matrices
# (described above).

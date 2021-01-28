from Bio.PDB import *
from pymol import cmd
import pandas as pd
import numpy as np
import pymol
import sys

def snps_to_aa(snps,imap_to_pdb):

	snps['tmp']=1
	imap_to_pdb['tmp']=1

	parser = PDBParser()
	structure_list = {}
	for pdb in set(imap_to_pdb['pdb'].values):
		pdbl=PDBList()
		pdbl.retrieve_pdb_file(pdb,pdir='ref/pdb',file_format='pdb',obsolete=False)
		structure = parser.get_structure(pdb,"ref/pdb/pdb%s.ent"%pdb)
		structure_list[pdb] = structure

	snps_pdb = pd.merge(snps, imap_to_pdb, on='tmp')

	for irow, row in snps_pdb.iterrows():
		entry = row['pdb']
		chain = row['chain']
		residue = int(row['res'])
		try:
			structure = structure_list.get(entry)
			atom = structure[0][chain][residue]["CA"]
			coord = atom.get_coord()
			snps_pdb.loc[irow,['x','y','z']] = coord
		except: continue

	snps_aa_complete = snps_pdb.dropna()
	uniq_map = snps_aa_complete.groupby(['pdb','chain'])['snp_id'].count().reset_index()
	iline = uniq_map['snp_id'].argmax()
	ipdb = uniq_map.loc[iline,'pdb']

	return snps_aa_complete[snps_aa_complete['pdb']==ipdb], ipdb

def plot(snps, pdb_id, active_sites, aa_sites, out, surface=False):

#	pymol.finish_launching()
	cmd.fetch(pdb_id)
	cmd.show_as('cartoon', pdb_id)
	cmd.color('white',pdb_id)

	for i, row in snps.iterrows():
		resi  = row['res']
		chain = row['chain']
		selec = 'snp%s'%i
		selec_atom = 'snp_atom%s'%i
		cmd.select(selec,'name ca and resi %s and chain %s'%(resi,chain))
		cmd.create(selec_atom, selec)
		cmd.set("sphere_scale",0.5)
		cmd.show('sphere',selec_atom)
		cmd.color('red',selec_atom)

	chains = np.unique(snps['chain'].values)
	# active state
	if active_sites:	
		for isite, site in enumerate(active_sites):
			for chain in chains:
				selec = 'active%s'%isite
				cmd.select(selec,'resi %s and chain %s'%(site,chain))
				if surface: 
					selec_surf = 'active_surf%s'%isite
					cmd.create(selec_surf,selec)
					cmd.show('surface',selec_surf)
					cmd.color('cyan', selec_surf)
					cmd.set('transparency', 0.3)
				else: cmd.color('cyan', selec)

	# aa state	
	if aa_sites:	
		for isite, site in enumerate(aa_sites):
			for chain in chains:	
				selec = 'aamod%s'%isite
				cmd.select(selec,'resi %s and chain %s'%(site,chain))
				if surface: 
					selec_surf = 'aamod_surf%s'%isite
					cmd.create(selec_surf,selec)
					cmd.show('surface',selec_surf)
					cmd.color('cyan', selec_surf)
					cmd.set('transparency', 0.3)
				else: cmd.color('cyan', selec)

	cmd.set('surface_quality',1)
	cmd.bg_color("white")
	cmd.zoom()
	cmd.png('%s.png'%out, dpi=300)
	cmd.save('%s.pse'%out)	

def main():

	gene = sys.argv[1]
	map_to_pdb_file = sys.argv[2]
	map_to_pdb = pd.read_csv(map_to_pdb_file,sep=' ',header=None,names=['uniprot','gene','pdb','chain'])
	snps_file = sys.argv[3]
	snps = pd.read_csv(snps_file,sep=' ',header=None,names=['snp_id','aa','res'])
	active_site = sys.argv[4]
	aa_site     = sys.argv[5]

	surface     = True if sys.argv[6] == 'surface' else False

	if active_site == 'None': active_site=None
	else: active_site = active_site.split(',')

	if aa_site == 'None': aa_site=None
	else: aa_site = aa_site.split(',')

	imap_to_pdb = map_to_pdb[map_to_pdb['gene']==gene]
	df, pdb = snps_to_aa(snps,imap_to_pdb)
	if surface: out = 'results/%s_%s_surf'%(gene,pdb)
	else: out = 'results/%s_%s'%(gene,pdb)
	plot(df,pdb,active_site, aa_site, out, surface=surface)


if __name__ == "__main__":
	main()

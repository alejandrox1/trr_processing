#ifndef file_pro_H
#define file_pro_H

void make_protein(std::vector<info>& atoms, std::vector<info*>& protein, unsigned int m, unsigned int n);

void mod_make_protein(std::vector<mod_info>& atoms, std::vector<mod_info*>& protein, unsigned int m, unsigned int n);

void ls_lh_gro(std::vector<std::string> &directory, unsigned int skip);

void parser(std::vector<info> &atoms, std::vector<mod_info> &mod_atoms, std::string file_handle, double &dimensions, int &num_atoms);

#endif

// Clipper Gemmi minimol test
// A test for MiniMol with Clipper-Gemmi interface

#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-gemmi.h> // for clipper::GemmiWriteOptions
#include <fstream>

using namespace clipper;

int main(int argc, char **argv)
{
    // Test clipper gemmi
    // test minimol
    MiniMol tempmol;
    
    GEMMIFile file;
    ::gemmi::Structure s; // for comparison
    try
    { // ground truth: reading structure using Gemmi function
        s = ::gemmi::read_structure_file(argv[1], ::gemmi::coor_format_from_ext(argv[1]));
        // casting to Clipper class
        GemmiStructure gs(s);
        std::cout << "Coordinate file : " << argv[1] << std::endl;
        std::cout << "Name : " << gs.name << std::endl;
        if (gs.spacegroup().is_null())
            std::cout << "  Spacegroup : NULL" << std::endl;
        else
            std::cout << "  Spacegroup : " << gs.spacegroup().symbol_hall() << std::endl;
        std::cout << "  Cell : " << gs.get_cell().format() << std::endl;

        switch (gs.input_format)
        {
        case ::gemmi::CoorFormat::Pdb:
            std::cout << "  Input format : PDB" << std::endl;
            break;
        case ::gemmi::CoorFormat::Mmcif:
            std::cout << "  Input format : Mmcif" << std::endl;
            break;
        case ::gemmi::CoorFormat::Unknown:
            std::cout << "  Input format : Unknown" << std::endl;
            break;
        }
        GemmiAtom_list temp_atmslist(gs.models[0].all(), 1);
        std::cout << "  Number of atoms in Structure : " << temp_atmslist.size() << std::endl;

        std::cout << "  Experimental method : " << gs.get_exptlmethod() << std::endl;
    }
    catch (std::runtime_error &e)
    {
        printf("File reading failed: %s\n", e.what());
        return 1;
    }
    try // test reading with GEMMIFile
    {
        // read file with GEMMIFile
        file.read_file(argv[1]);
        // import minimol
        file.import_minimol(tempmol);
        Atom_list atoms = tempmol.atom_list();
        std::cout << "Coordinate file : " << argv[1] << std::endl;
        std::cout << "  Spacegroup : " << file.spacegroup().symbol_hall() << std::endl;
        std::cout << "  Cell : " << file.cell().format() << std::endl;
        std::cout << "  Number of atoms read : " << atoms.size() << std::endl;

        // get GemmiStructure
        GemmiStructure temp_st;
        temp_st = file.get_gemmi_structure();
        std::cout << "Structure from GEMMIFile : " << std::endl;
        std::cout << "  Name : " << temp_st.name << std::endl;
        GemmiAtom_list temp_atmslist(temp_st.models[0].all(), 1);
        std::cout << "  Number of atoms in Structure : " << temp_atmslist.size() << std::endl;

        // Mat33<float> m;
        Vec3<> v(2.0, 2.0, 2.0);
        RTop_orth rt(Mat33<>::identity(), v);
        MiniMol mol;
        mol.model() = tempmol.model();
        for (int c = 0; c < mol.size(); c++)
            for (int r = 0; r < mol[c].size(); r++)
                for (int a = 0; a < mol[c][r].size(); a++)
                {
                    std::cout << mol[c].id() << "\t" << mol[c][r].id() << "\t" << mol[c][r].type() << "\t" << mol[c][r][a].id() << ">";
                    std::cout << "\t" << dynamic_cast<const Property<String> &>(mol[c][r][a].get_property("CID")).value();
                    // std::cout << "\t>> " << dynamic_cast<const Property<String> &>(mol2[p][m][a].get_property("CID")).value();
                    if (mol[c][r][a].exists_property("AltConf"))
                    {

                        std::cout << "\t" << dynamic_cast<const Property<String> &>(mol[c][r][a].get_property("AltConf")).value();
                    }
                    std::cout << " ANISO : " << bool(mol[c][r][a].u_aniso_orth().is_null());
                    std::cout << "\n";
                    // try translating coordinates with transform
                    mol[c][r][a].transform(rt);
                }
        // write files
        // export minimol
        std::cout << "Writing coordinate files " << std::endl;
        file.export_minimol(mol);
        GemmiWriteOptions write_options;
        file.write_file("cgtest_mod.pdb"); // write default pdb
        // write pdb with minimal data
        write_options.Minimal = true;
        file.write_file("cgtest_mod_minimal.pdb", GEMMIFile::PDB, write_options);
        // write Pdbx style cif with minimal data
        file.set_cif_style(GEMMIFile::Pdbx);
        file.write_file("cgtest_mod_minimal_pdbx.cif", GEMMIFile::CIF, write_options);
        // write aligned cif with minimal data
        file.set_cif_style(GEMMIFile::Aligned);
        file.write_file("cgtest_mod_minimal_aligned.cif", GEMMIFile::CIF, write_options);

        // update cif doc, write files based on groups set
        write_options.Minimal = false;
        write_options.UpdateCifDoc = true;
        GemmiMmcifOutputGroups groups(true);
        file.set_mmcif_output_groups(groups);
        file.write_file("cgtest_groupstrue.cif", GEMMIFile::CIF, write_options);

        // test exporting minimol into new GEMMIFile object and write file
        GEMMIFile file2;
        file2.export_minimol(mol);
        file2.write_file("cgtest_cpy.pdb");
        file2.set_mmcif_output_groups(groups);
        file2.write_file("cgtest_cpy.cif", GEMMIFile::CIF, write_options);
    }
    catch (Message_fatal)
    {
        std::cout << "FAILED TO READ COORDINATE FILE " << argv[1] << '\n';
        return 1;
    }
}

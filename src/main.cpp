#include <iostream>
#include <string>
#include "gmshparsercpp.h"
#include "exodusIIcpp.h"

/// Mesh dimension
int dim = -1;
int side_set_dim = -1;
int node_set_dim = -1;

std::vector<double> x;
std::vector<double> y;
std::vector<double> z;

std::map<int, int> node_map;

std::vector<exodusIIcpp::ElementBlock> element_blocks;
std::vector<exodusIIcpp::SideSet> side_sets;
std::vector<exodusIIcpp::NodeSet> node_sets;

// NOTE: this may be needed per block
std::map<int, int> elem_map;

/// GMSH elem type to exodusII string type representation
std::map<int, const char *> exo_elem_type = {
        {1, "EDGE2"},
        {2, "TRI3"},
        {3, "QUAD4"},
        {4, "TET4" },
        {5, "HEX8" }
};

/// GMSH elem type to number of nodes per that element type
std::map<int, int> nodes_per_elem {
        {1, 2},
        {2, 3},
        {3, 4},
        {4, 4},
        {5, 8}
};

//

void print_help() {
    std::cout << "gmsh2exo <input> <output>" << std::endl;
}

void read_gmsh_file(const std::string &file_name) {
    gmshparsercpp::MshFile f(file_name);
    f.parse();

    const std::vector<gmshparsercpp::MshFile::PhysicalName> &phys_entities = f.get_physical_names();
    std::vector<std::vector<const gmshparsercpp::MshFile::PhysicalName *>> phys_ents_dim;
    phys_ents_dim.resize(4);
    for (const auto &pe: phys_entities)
        phys_ents_dim[pe.dimension].push_back(&pe);

    const std::vector<gmshparsercpp::MshFile::ElementBlock> &el_blks = f.get_element_blocks();
    std::vector<std::vector<const gmshparsercpp::MshFile::ElementBlock *>> el_blk_dim;
    el_blk_dim.resize(4);
    for (const auto &eb: el_blks)
        el_blk_dim[eb.dimension].push_back(&eb);

    // the element block with the highest dimension that has something in it will be the final mesh dimension
    for (int i = 0; i < 4; i++)
        if (el_blk_dim[i].size() > 0)
            dim = i;
    side_set_dim = dim - 1;
    node_set_dim = 0;

    const std::vector<gmshparsercpp::MshFile::Node> &nodes = f.get_nodes();
//    std::cout << "nodes = " << nodes.size() << std::endl;
    for (const auto &nd: nodes) {
//        std::cout << " - dim = " << nd.dimension;
//        std::cout << ", ent_tag = " << nd.entity_tag;
//        std::cout << ", tags = ";
//        for (const auto &tag: nd.tags)
//            std::cout << " " << tag;
//        std::cout << ", coords = ";
//        for (const auto &c: nd.coordinates)
//            std::cout << " (" << c.x << ", " << c.y << ", " << c.z << ")";
//        std::cout << std::endl;

        if (!nd.tags.empty()) {
            for (int j = 0; j < nd.tags.size(); j++) {
                int local_id = node_map.size();
                const auto &id = nd.tags[j];
                node_map[local_id] = id;

                const auto &c = nd.coordinates[j];
                x.push_back(c.x);
                y.push_back(c.y);
                z.push_back(c.z);
            }
        }
    }

//    std::cout << "elems" << std::endl;
    for (const auto &eb: el_blks) {
//        std::cout << "- dim = " << eb.dimension;
//        std::cout << ", el_type = " << eb.element_type;
//        std::cout << ", tag = " << eb.tag;
//        std::cout << std::endl;
//        std::cout << "  elems:" << std::endl;
//        for (const auto &e: eb.elements) {
//            std::cout << "  - tag = " << e.tag;
//            std::cout << ", nodes =";
//            for (const auto &nt: e.node_tags)
//                std::cout << " " << nt;
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;

        if (eb.dimension == dim) {
            exodusIIcpp::ElementBlock exo_eb;
            exo_eb.set_id(eb.tag);
            std::vector<int> connect;
            for (const auto &elem: eb.elements)
                for (const auto & nid: elem.node_tags)
                    connect.push_back(nid);
            int n_nodes_pre_elem = nodes_per_elem[eb.element_type];
            exo_eb.set_connectivity(exo_elem_type.at(eb.element_type), eb.elements.size(), n_nodes_pre_elem, connect);
            element_blocks.push_back(exo_eb);
        }
        else if (eb.dimension == side_set_dim) {
            // exodusIIcpp::SideSet ss;
            // // FIXME: this should be the tag from the physical entity corresponding to this element block
            // ss.set_id(side_sets.size() + 1);
            // // FIXME: this should be the physical entity name corresponding to this element block
            // ss.set_name();
            // std::vector<int> elem_ids;
            // std::vector<int> side_ids;
            //
            // side_sets.push_back(ss);
        }
        else if (eb.dimension == node_set_dim) {
            exodusIIcpp::NodeSet ns;
            ns.set_id(node_sets.size() + 1);
            // TODO: this should be set to the physical entity name corresponding to this element block
            // ns.set_name();
            std::vector<int> nodes;
            for (const auto &elem: eb.elements)
                for (const auto & nid: elem.node_tags)
                    nodes.push_back(nid);
            ns.set_nodes(nodes);
            node_sets.push_back(ns);
        }
    }
}

void write_exodus_file(const std::string &file_name) {
    exodusIIcpp::File f(file_name, exodusIIcpp::FileAccess::WRITE);
    int n_nodes = x.size();
    int n_elems = 0;
    for (const auto & eb: element_blocks)
        n_elems += eb.get_num_elements();
    int n_elem_blks = element_blocks.size();
    int n_node_sets = 0;
    int n_side_sets = 0;

    f.init("", dim, n_nodes, n_elems, n_elem_blks, n_node_sets, n_side_sets);
//    std::cout << "dim = " << dim << ", n_elems = " << n_elems << ", n_elem_blks = " << n_elem_blks << std::endl;
    if (dim == 1)
        f.write_coords(x);
    else if (dim == 2)
        f.write_coords(x, y);
    else
        f.write_coords(x, y, z);
    f.write_coord_names();

    for (const auto &eb: element_blocks)
        f.write_block(eb.get_id(), eb.get_element_type().c_str(), eb.get_num_elements(), eb.get_connectivity());

    for (const auto & ss: side_sets)
        f.write_side_set(ss.get_id(), ss.get_element_ids(), ss.get_side_ids());
//    f.write_side_set_names();

    for (const auto & ns: node_sets)
        f.write_node_set(ns.get_id(), ns.get_node_ids());
//    f.write_node_set_names();

    f.close();
}

void convert(const std::string &input_file_name, const std::string &output_file_name) {
    read_gmsh_file(input_file_name);
    write_exodus_file(output_file_name);
}

int main(int argc, char *argv[]) {
    try {
        if (argc > 2) {
            convert(argv[1], argv[2]);
        } else
            print_help();
    }
    catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

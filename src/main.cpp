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
// (node_ids) -> (elem, side)
std::map<std::vector<int>, std::pair<int, int>> elem_sides;
std::map<int, const gmshparsercpp::MshFile::PhysicalName *> phys_ent_by_tag;

std::vector<exodusIIcpp::ElementBlock> element_blocks;
std::map<int, exodusIIcpp::SideSet> side_sets;
std::map<int, exodusIIcpp::NodeSet> node_sets;

// NOTE: this may be needed per block
std::map<int, int> elem_map;

/// GMSH elem type to exodusII string type representation
std::map<int, const char *> exo_elem_type = {
        {1, "EDGE2"},
        {2, "TRI3"},
        {3, "QUAD4"},
        {4, "TET4"},
        {5, "HEX8"}
};

/// GMSH elem type to number of nodes per that element type
std::map<int, int> nodes_per_elem{
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

std::vector<int> build_side_key_edge2(int id1, int id2) {
    std::vector<int> key;
    if (id1 < id2) {
        key.push_back(id1);
        key.push_back(id2);
    } else {
        key.push_back(id2);
        key.push_back(id1);
    }
    return key;
}

void read_physical_entities(const std::vector<gmshparsercpp::MshFile::PhysicalName> &phys_entities) {
    for (const auto &pe: phys_entities) {
        phys_ent_by_tag[pe.tag] = &pe;
    }
}

void analyze_mesh(const std::vector<gmshparsercpp::MshFile::ElementBlock> &el_blks) {
    std::vector<std::vector<const gmshparsercpp::MshFile::ElementBlock *>> el_blk_dim;
    el_blk_dim.resize(4);
    for (const auto &eb: el_blks)
        el_blk_dim[eb.dimension].push_back(&eb);

    // the element block with the highest dimension that has something in it will be the final mesh dimension
    for (int i = 0; i < 4; i++)
        if (!el_blk_dim[i].empty())
            dim = i;
    side_set_dim = dim - 1;
    node_set_dim = 0;
}

void build_coordinates(const std::vector<gmshparsercpp::MshFile::Node> &nodes) {
    for (const auto &nd: nodes) {
        if (!nd.tags.empty()) {
            for (int j = 0; j < nd.tags.size(); j++) {
                int local_id = (int) node_map.size();
                const auto &id = nd.tags[j];
                node_map[local_id] = id;

                const auto &c = nd.coordinates[j];
                x.push_back(c.x);
                y.push_back(c.y);
                z.push_back(c.z);
            }
        }
    }
}

void build_element_blocks(const std::vector<gmshparsercpp::MshFile::ElementBlock> &el_blks,
                          const std::vector<gmshparsercpp::MshFile::MultiDEntity> &entities) {
    std::map<int, const gmshparsercpp::MshFile::MultiDEntity *> ents_by_id;
    for (const auto &ent: entities)
        ents_by_id[ent.tag] = &ent;

    unsigned int eid = 0;
    for (const auto &eb: el_blks) {
        if (eb.dimension == dim) {
            exodusIIcpp::ElementBlock exo_eb;

            const auto &ent = ents_by_id[eb.tag];
            if (!ent->physical_tags.empty()) {
                auto id = ent->physical_tags[0];
                exo_eb.set_name(phys_ent_by_tag[id]->name);
                exo_eb.set_id(id);
            } else
                exo_eb.set_id(eb.tag);

            std::vector<int> connect;
            for (const auto &elem: eb.elements) {
                for (const auto &nid: elem.node_tags)
                    connect.push_back(nid);

                auto n_node_tags = elem.node_tags.size();
                for (unsigned int i = 0; i < n_node_tags; i++) {
                    std::vector<int> side_key;
                    if (dim == 1) {
                        side_key.push_back(elem.node_tags[i]);
                    } else if (dim == 2) {
                        side_key = build_side_key_edge2(elem.node_tags[i], elem.node_tags[(i + 1) % n_node_tags]);
                    } else if (dim == 3) {
                        throw std::runtime_error("not implemented yet");
                    }
                    elem_sides[side_key] = std::pair(eid + 1, i + 1);
                }
                eid++;
            }
            int n_nodes_pre_elem = nodes_per_elem[eb.element_type];
            exo_eb.set_connectivity(exo_elem_type.at(eb.element_type), (int) eb.elements.size(), n_nodes_pre_elem, connect);
            element_blocks.push_back(exo_eb);
        }
    }
}

void build_side_sets(const std::vector<gmshparsercpp::MshFile::ElementBlock> &el_blks,
                     const std::vector<gmshparsercpp::MshFile::MultiDEntity> &entities) {
    std::map<int, const gmshparsercpp::MshFile::MultiDEntity *> ents_by_id;
    for (const auto &ent: entities) {
        if (!ent.physical_tags.empty()) {
            for (const auto &id: ent.physical_tags) {
                exodusIIcpp::SideSet ss;
                ss.set_id(id);
                ss.set_name(phys_ent_by_tag[id]->name);
                side_sets[id] = ss;
            }
        }
        ents_by_id[ent.tag] = &ent;
    }

    for (const auto &eb: el_blks) {
        if (eb.dimension == side_set_dim) {
            const auto &ent = ents_by_id[eb.tag];

            for (const auto &id: ent->physical_tags) {
                auto &ss = side_sets[id];

                if (eb.element_type == 1) {
                    // 2-node edge
                    for (const auto &elem: eb.elements) {
                        std::vector<int> side_key = build_side_key_edge2(elem.node_tags[0], elem.node_tags[1]);
                        const auto el_side_pair = elem_sides[side_key];
                        ss.add(el_side_pair.first, el_side_pair.second);
                    }
                }
            }
        }
    }
}

void read_gmsh_file(const std::string &file_name) {
    gmshparsercpp::MshFile f(file_name);
    f.parse();

    const std::vector<gmshparsercpp::MshFile::PhysicalName> &phys_entities = f.get_physical_names();
    read_physical_entities(phys_entities);

    const std::vector<gmshparsercpp::MshFile::Node> &nodes = f.get_nodes();
    const std::vector<gmshparsercpp::MshFile::ElementBlock> &el_blks = f.get_element_blocks();
    analyze_mesh(el_blks);
    build_coordinates(nodes);
    if (dim == 2) {
        const std::vector<gmshparsercpp::MshFile::MultiDEntity> &ents = f.get_surface_entities();
        build_element_blocks(el_blks, ents);
    }

    if (side_set_dim == 1) {
        const std::vector<gmshparsercpp::MshFile::MultiDEntity> &curve_ents = f.get_curve_entities();
        build_side_sets(el_blks, curve_ents);
    } else {
        throw std::runtime_error("side sets with dim other than 1 are not implemented yet");
    }
}

void write_exodus_coordinates(exodusIIcpp::File &f) {
    if (dim == 1)
        f.write_coords(x);
    else if (dim == 2)
        f.write_coords(x, y);
    else
        f.write_coords(x, y, z);
    f.write_coord_names();
}

void write_exodus_element_blocks(exodusIIcpp::File &f) {
    std::vector<std::string> block_names;
    for (const auto &eb: element_blocks) {
        f.write_block(eb.get_id(), eb.get_element_type().c_str(), eb.get_num_elements(), eb.get_connectivity());
        block_names.push_back(eb.get_name());
    }
    f.write_block_names(block_names);
}

void write_exodus_side_sets(exodusIIcpp::File &f) {
    std::vector<std::string> sideset_names;
    for (const auto &it: side_sets) {
        const auto &ss = it.second;
        if (ss.get_size() > 0) {
            f.write_side_set(ss.get_id(), ss.get_element_ids(), ss.get_side_ids());
            sideset_names.push_back(ss.get_name());
        }
    }
    if (!sideset_names.empty())
        f.write_side_set_names(sideset_names);
}

void write_exodus_node_sets(exodusIIcpp::File &f) {
    std::vector<std::string> nodeset_names;
    for (const auto &it: node_sets) {
        const auto &ns = it.second;
        if (ns.get_size() > 0) {
            f.write_node_set(ns.get_id(), ns.get_node_ids());
            nodeset_names.push_back(ns.get_name());
        }
    }
    if (!nodeset_names.empty())
        f.write_node_set_names(nodeset_names);
}

void write_exodus_file(const std::string &file_name) {
    exodusIIcpp::File f(file_name, exodusIIcpp::FileAccess::WRITE);

    int n_nodes = (int) x.size();
    int n_elems = 0;
    for (const auto &eb: element_blocks)
        n_elems += eb.get_num_elements();
    int n_elem_blks = (int) element_blocks.size();
    int n_node_sets = (int) node_sets.size();
    int n_side_sets = (int) side_sets.size();
    f.init("", dim, n_nodes, n_elems, n_elem_blks, n_node_sets, n_side_sets);

    write_exodus_coordinates(f);
    write_exodus_element_blocks(f);
    write_exodus_side_sets(f);
    write_exodus_node_sets(f);

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

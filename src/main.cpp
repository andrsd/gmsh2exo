#include <iostream>
#include <string>
#include "gmshparsercpp.h"
#include "exodusIIcpp.h"
#include "cxxopts.hpp"
#include "Gmsh2ExoConfig.h"

/// Mesh dimension
int dim = -1;
int side_set_dim = -1;
int node_set_dim = -1;

std::vector<double> x;
std::vector<double> y;
std::vector<double> z;

std::map<int, int> node_map;
/// Blocks per dimension
std::vector<std::vector<const gmshparsercpp::MshFile::ElementBlock *>> el_blk_dim;
// (node_ids) -> (elem, side)
std::map<std::vector<int>, std::pair<int, int>> elem_sides;
std::map<int, const gmshparsercpp::MshFile::PhysicalName *> phys_ent_by_tag;

std::vector<exodusIIcpp::ElementBlock> element_blocks;
std::map<int, exodusIIcpp::SideSet> side_sets;
std::map<int, exodusIIcpp::NodeSet> node_sets;

// NOTE: this may be needed per block
std::map<int, int> elem_map;

namespace gmsh {
enum ElementType { EDGE2 = 1, TRI3 = 2, QUAD4 = 3, TET4 = 4, HEX8 = 5 };
}

/// GMSH elem type to exodusII string type representation
std::map<int, const char *> exo_elem_type = { { gmsh::EDGE2, "EDGE2" },
                                              { gmsh::TRI3, "TRI3" },
                                              { gmsh::QUAD4, "QUAD4" },
                                              { gmsh::TET4, "TET4" },
                                              { gmsh::HEX8, "HEX8" } };

/// GMSH elem type to number of nodes per that element type
std::map<int, int> nodes_per_elem = { { gmsh::EDGE2, 2 },
                                      { gmsh::TRI3, 3 },
                                      { gmsh::QUAD4, 4 },
                                      { gmsh::TET4, 4 },
                                      { gmsh::HEX8, 8 } };

// node ordering from GMSH to exodusII
std::vector<std::vector<int>> node_order = {
    {}, { 0, 1 }, { 0, 1, 2 }, { 0, 1, 2, 3 }, { 2, 1, 0, 3 }, { 0, 1, 2, 3, 4, 5, 6, 7 }
};

//

std::vector<int>
build_side_key_edge2(int id1, int id2)
{
    std::vector<int> key;
    if (id1 < id2) {
        key.push_back(id1);
        key.push_back(id2);
    }
    else {
        key.push_back(id2);
        key.push_back(id1);
    }
    return key;
}

std::vector<int>
build_side_key_tri3(int id1, int id2, int id3)
{
    std::vector<int> key = { id1, id2, id3 };
    std::sort(key.begin(), key.end());
    return key;
}

std::vector<int>
build_side_key_quad4(int id1, int id2, int id3, int id4)
{
    std::vector<int> key = { id1, id2, id3, id4 };
    std::sort(key.begin(), key.end());
    return key;
}

void
read_physical_entities(const std::vector<gmshparsercpp::MshFile::PhysicalName> & phys_entities)
{
    for (const auto & pe : phys_entities) {
        phys_ent_by_tag[pe.tag] = &pe;
    }
}

void
analyze_mesh(const std::vector<gmshparsercpp::MshFile::ElementBlock> & el_blks)
{
    el_blk_dim.resize(4);
    for (const auto & eb : el_blks)
        el_blk_dim[eb.dimension].push_back(&eb);

    // the element block with the highest dimension that has something in it will be the final mesh
    // dimension
    for (int i = 0; i < 4; i++)
        if (!el_blk_dim[i].empty())
            dim = i;
    side_set_dim = dim - 1;
    node_set_dim = 0;
}

void
build_coordinates(const std::vector<gmshparsercpp::MshFile::Node> & nodes)
{
    for (const auto & nd : nodes) {
        if (!nd.tags.empty()) {
            for (int j = 0; j < nd.tags.size(); j++) {
                int local_id = (int) node_map.size();
                const auto & id = nd.tags[j];
                node_map[local_id] = id;

                const auto & c = nd.coordinates[j];
                x.push_back(c.x);
                y.push_back(c.y);
                z.push_back(c.z);
            }
        }
    }
}

void
build_element_blocks(const std::vector<const gmshparsercpp::MshFile::ElementBlock *> & el_blks,
                     const std::vector<gmshparsercpp::MshFile::MultiDEntity> & entities)
{
    std::map<int, const gmshparsercpp::MshFile::MultiDEntity *> ents_by_id;
    for (const auto & ent : entities)
        ents_by_id[ent.tag] = &ent;

    unsigned int eid = 0;
    for (const auto & eb : el_blks) {
        exodusIIcpp::ElementBlock exo_eb;

        const auto & ent = ents_by_id[eb->tag];
        if (!ent->physical_tags.empty()) {
            auto id = ent->physical_tags[0];
            auto it = phys_ent_by_tag.find(id);
            if (it != phys_ent_by_tag.end())
                exo_eb.set_name(it->second->name);
            exo_eb.set_id(id);
        }
        else
            exo_eb.set_id(eb->tag);

        std::vector<int> connect;
        for (const auto & elem : eb->elements) {
            std::vector<int> el_nodes;
            auto n_node_tags = elem.node_tags.size();
            for (std::size_t i = 0; i < n_node_tags; i++) {
                const auto & nid = elem.node_tags[node_order[eb->element_type][i]];
                el_nodes.push_back(nid);
            }
            connect.insert(connect.end(), el_nodes.begin(), el_nodes.end());

            if (dim == 1) {
                for (unsigned int i = 0; i < n_node_tags; i++) {
                    std::vector<int> side_key = { el_nodes[i] };
                    elem_sides[side_key] = std::pair(eid + 1, i + 1);
                }
            }
            else if (dim == 2) {
                for (unsigned int i = 0; i < n_node_tags; i++) {
                    std::vector<int> side_key =
                        build_side_key_edge2(el_nodes[i], el_nodes[(i + 1) % n_node_tags]);
                    elem_sides[side_key] = std::pair(eid + 1, i + 1);
                }
            }
            else if (dim == 3) {
                // Note: see Sjaardema, G. D., Schoof, L. A. & Yarberry, V. R. EXODUS: A Finite
                // Element Data Model. 148 (2019) for how sides are numbered on different elements
                // (fig 4.15, pp. 28)
                if (eb->element_type == gmsh::TET4) {
                    std::vector<std::vector<int>> sides = { { 0, 1, 3 },
                                                            { 1, 2, 3 },
                                                            { 0, 2, 3 },
                                                            { 0, 1, 2 } };
                    for (std::size_t s = 0; s < sides.size(); s++) {
                        std::vector<int> side_key;
                        side_key = build_side_key_tri3(el_nodes[sides[s][0]],
                                                       el_nodes[sides[s][1]],
                                                       el_nodes[sides[s][2]]);
                        elem_sides[side_key] = std::pair(eid + 1, s + 1);
                    }
                }
                else if (eb->element_type == gmsh::HEX8) {
                    std::vector<std::vector<int>> sides = { { 0, 1, 5, 4 }, { 1, 2, 6, 5 },
                                                            { 3, 2, 6, 7 }, { 0, 3, 7, 4 },
                                                            { 0, 1, 2, 3 }, { 4, 5, 6, 7 } };
                    for (std::size_t s = 0; s < sides.size(); s++) {
                        std::vector<int> side_key;
                        side_key = build_side_key_quad4(el_nodes[sides[s][0]],
                                                        el_nodes[sides[s][1]],
                                                        el_nodes[sides[s][2]],
                                                        el_nodes[sides[s][3]]);
                        elem_sides[side_key] = std::pair(eid + 1, s + 1);
                    }
                }
                else
                    throw std::runtime_error("not implemented yet");
            }
            eid++;
        }
        int n_nodes_pre_elem = nodes_per_elem[eb->element_type];
        exo_eb.set_connectivity(exo_elem_type.at(eb->element_type),
                                (int) eb->elements.size(),
                                n_nodes_pre_elem,
                                connect);
        element_blocks.push_back(exo_eb);
    }
}

void
build_side_sets(const std::vector<const gmshparsercpp::MshFile::ElementBlock *> & el_blks,
                const std::vector<gmshparsercpp::MshFile::MultiDEntity> & entities)
{
    std::map<int, const gmshparsercpp::MshFile::MultiDEntity *> ents_by_id;
    for (const auto & ent : entities) {
        if (!ent.physical_tags.empty()) {
            for (const auto & id : ent.physical_tags) {
                exodusIIcpp::SideSet ss;
                ss.set_id(id);
                auto it = phys_ent_by_tag.find(id);
                if (it != phys_ent_by_tag.end())
                    ss.set_name(it->second->name);
                side_sets[id] = ss;
            }
        }
        ents_by_id[ent.tag] = &ent;
    }

    for (const auto & eb : el_blks) {
        const auto & ent = ents_by_id[eb->tag];

        for (const auto & id : ent->physical_tags) {
            auto & ss = side_sets[id];

            for (const auto & elem : eb->elements) {
                std::vector<int> side_key;
                switch (eb->element_type) {
                case gmsh::EDGE2:
                    side_key = build_side_key_edge2(elem.node_tags[0], elem.node_tags[1]);
                    break;
                case gmsh::TRI3:
                    side_key = build_side_key_tri3(elem.node_tags[0],
                                                   elem.node_tags[1],
                                                   elem.node_tags[2]);
                    break;
                case gmsh::QUAD4:
                    side_key = build_side_key_quad4(elem.node_tags[0],
                                                    elem.node_tags[1],
                                                    elem.node_tags[2],
                                                    elem.node_tags[3]);
                    break;
                }
                const auto el_side_pair = elem_sides[side_key];
                ss.add(el_side_pair.first, el_side_pair.second);
            }
        }
    }
}

void
read_gmsh_file(const std::string & file_name)
{
    gmshparsercpp::MshFile f(file_name);
    f.parse();

    const std::vector<gmshparsercpp::MshFile::PhysicalName> & phys_entities =
        f.get_physical_names();
    read_physical_entities(phys_entities);

    const std::vector<gmshparsercpp::MshFile::Node> & nodes = f.get_nodes();
    const std::vector<gmshparsercpp::MshFile::ElementBlock> & el_blks = f.get_element_blocks();
    analyze_mesh(el_blks);
    build_coordinates(nodes);

    const std::vector<gmshparsercpp::MshFile::MultiDEntity> * ents = nullptr;
    switch (dim) {
    case 1:
        ents = &f.get_curve_entities();
        break;
    case 2:
        ents = &f.get_surface_entities();
        break;
    case 3:
        ents = &f.get_volume_entities();
        break;
    default:
        throw std::runtime_error("Unsupported dimension.");
    }
    build_element_blocks(el_blk_dim[dim], *ents);

    const std::vector<gmshparsercpp::MshFile::MultiDEntity> * sideset_ents = nullptr;
    switch (side_set_dim) {
    case 1:
        sideset_ents = &f.get_curve_entities();
        break;
    case 2:
        sideset_ents = &f.get_surface_entities();
        break;
    default:
        throw std::runtime_error("Unsupported side sets dim.");
    }
    build_side_sets(el_blk_dim[side_set_dim], *sideset_ents);
}

void
write_exodus_coordinates(exodusIIcpp::File & f)
{
    if (dim == 1)
        f.write_coords(x);
    else if (dim == 2)
        f.write_coords(x, y);
    else
        f.write_coords(x, y, z);
    f.write_coord_names();
}

void
write_exodus_element_blocks(exodusIIcpp::File & f)
{
    std::vector<std::string> block_names;
    for (const auto & eb : element_blocks) {
        f.write_block(eb.get_id(),
                      eb.get_element_type().c_str(),
                      eb.get_num_elements(),
                      eb.get_connectivity());
        block_names.push_back(eb.get_name());
    }
    f.write_block_names(block_names);
}

void
write_exodus_side_sets(exodusIIcpp::File & f)
{
    std::vector<std::string> sideset_names;
    for (const auto & it : side_sets) {
        const auto & ss = it.second;
        if (ss.get_size() > 0) {
            f.write_side_set(ss.get_id(), ss.get_element_ids(), ss.get_side_ids());
            sideset_names.push_back(ss.get_name());
        }
    }
    if (!sideset_names.empty())
        f.write_side_set_names(sideset_names);
}

void
write_exodus_node_sets(exodusIIcpp::File & f)
{
    std::vector<std::string> nodeset_names;
    for (const auto & it : node_sets) {
        const auto & ns = it.second;
        if (ns.get_size() > 0) {
            f.write_node_set(ns.get_id(), ns.get_node_ids());
            nodeset_names.push_back(ns.get_name());
        }
    }
    if (!nodeset_names.empty())
        f.write_node_set_names(nodeset_names);
}

void
write_exodus_file(const std::string & file_name)
{
    exodusIIcpp::File f(file_name, exodusIIcpp::FileAccess::WRITE);

    int n_nodes = (int) x.size();
    int n_elems = 0;
    for (const auto & eb : element_blocks)
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

void
convert(const std::string & input_file_name, const std::string & output_file_name)
{
    read_gmsh_file(input_file_name);
    write_exodus_file(output_file_name);
}

int
main(int argc, char * argv[])
{
    cxxopts::Options options("gmsh2exo", "Convert GMSH mesh files into exodusII files");

    // clang-format off
    options.add_options()
        ("help", "Show this help page")
        ("v,version", "Show the version")
        ("gmsh_file", "GMSH file name", cxxopts::value<std::string>())
        ("exo2_file", "ExodusII file name", cxxopts::value<std::string>())
    ;
    options.parse_positional({ "gmsh_file", "exo2_file" });
    options.positional_help("<gmsh file> <exodusII file>");
    // clang-format on

    cxxopts::ParseResult result;
    try {
        result = options.parse(argc, argv);
    }
    catch (const cxxopts::exceptions::exception & e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cout << options.help();
        return 1;
    }

    if (result.count("version"))
        std::cout << "gmsh2exo version " << GMSH2EXO_VERSION << std::endl;
    else if (result.count("help"))
        std::cout << options.help();
    else if (result.count("gmsh_file") && result.count("exo2_file")) {
        try {
            convert(result["gmsh_file"].as<std::string>(), result["exo2_file"].as<std::string>());
        }
        catch (std::exception & e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
    }
    else
        std::cout << options.help();

    return 0;
}

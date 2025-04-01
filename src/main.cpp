// SPDX-FileCopyrightText: 2023 (c) David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

// References:
// [1] Sjaardema, G. D., Schoof, L. A. & Yarberry, V. R. EXODUS: A Finite Element Data Model. 148
// (2019)

#include "gmshparsercpp/Enums.h"
#include "gmshparsercpp/MshFile.h"
#include "exodusIIcpp/exodusIIcpp.h"
#include "cxxopts.hpp"
#include "fmt/format.h"
#include <string>
#include <array>

/// Mesh dimension
int dim = -1;
/// Side set dimension
int side_set_dim = -1;
/// Node set dimension
int node_set_dim = -1;

/// x-coordinates
std::vector<double> x;
/// y-coordinates
std::vector<double> y;
/// z-coordinates
std::vector<double> z;
// Map from GMSH node id to node id (i.e. index into coordinates, but not 1-based exodusII index)
std::map<int, int> node_map;
/// Blocks per dimension
std::vector<std::vector<const gmshparsercpp::MshFile::ElementBlock *>> el_blk_dim;
// (node_ids) -> (elem, side)
std::map<std::vector<int>, std::pair<int, int>> elem_sides;
/// Map from physical entity tag to physical entity
std::map<int, const gmshparsercpp::MshFile::PhysicalName *> phys_ent_by_tag;

/// Element blocks to be written to exodusII file
std::vector<exodusIIcpp::ElementBlock> element_blocks;
/// Side sets to be written to exodusII file
std::map<int, exodusIIcpp::SideSet> side_sets;
/// Node sets to be written to exodusII file
std::map<int, exodusIIcpp::NodeSet> node_sets;

// NOTE: this may be needed per block
std::map<int, int> elem_map;

/// Bounding box
struct BoundingBox {
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;

    BoundingBox() :
        xmin(std::numeric_limits<double>::max()),
        xmax(std::numeric_limits<double>::lowest()),
        ymin(std::numeric_limits<double>::max()),
        ymax(std::numeric_limits<double>::lowest()),
        zmin(std::numeric_limits<double>::max()),
        zmax(std::numeric_limits<double>::lowest())
    {
    }
};

/// GMSH elem type to exodusII string type representation
std::map<int, const char *> exo_elem_type = { { gmshparsercpp::LINE2, "EDGE2" },
                                              { gmshparsercpp::TRI3, "TRI3" },
                                              { gmshparsercpp::QUAD4, "QUAD4" },
                                              { gmshparsercpp::TET4, "TET4" },
                                              { gmshparsercpp::HEX8, "HEX8" } };

/// GMSH elem type to number of nodes per that element type
std::map<int, int> nodes_per_elem = { { gmshparsercpp::LINE2, 2 },
                                      { gmshparsercpp::TRI3, 3 },
                                      { gmshparsercpp::QUAD4, 4 },
                                      { gmshparsercpp::TET4, 4 },
                                      { gmshparsercpp::HEX8, 8 } };

// node ordering from GMSH to exodusII
std::vector<std::vector<int>> node_order = {
    {}, { 0, 1 }, { 0, 1, 2 }, { 0, 1, 2, 3 }, { 2, 1, 0, 3 }, { 0, 1, 2, 3, 4, 5, 6, 7 }
};

//

/// Build a side set key for EDGE2 elements
///
/// @param id1 First node id
/// @param id2 Second node id
/// @return Side set key
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

/// Build a side set key for a TRI3 element
///
/// @param id1 First node id
/// @param id2 Second node id
/// @param id3 Third node id
/// @return Side set key
std::vector<int>
build_side_key_tri3(int id1, int id2, int id3)
{
    std::vector<int> key = { id1, id2, id3 };
    std::sort(key.begin(), key.end());
    return key;
}

/// Build a side set key for a QUAD4 element
///
/// @param id1 First node id
/// @param id2 Second node id
/// @param id3 Third node id
/// @param id4 Fourth node id
/// @return Side set key
std::vector<int>
build_side_key_quad4(int id1, int id2, int id3, int id4)
{
    std::vector<int> key = { id1, id2, id3, id4 };
    std::sort(key.begin(), key.end());
    return key;
}

/// Build a physical entity to tag map
///
/// This function builds a map from physical entity tag to physical entity
///
/// @param phys_entities Physical entities
void
read_physical_entities(const std::vector<gmshparsercpp::MshFile::PhysicalName> & phys_entities)
{
    for (const auto & pe : phys_entities) {
        phys_ent_by_tag[pe.tag] = &pe;
    }
}

/// Compute mesh bounding box
///
/// This function computes the bounding box of the mesh
///
/// @return Mesh bounding box
BoundingBox
compute_mesh_bounding_box()
{
    // this assumes that `x`, `y` and `z` contain valid coordinate data
    BoundingBox bbox;
    if ((x.size() == y.size()) and (y.size() == z.size())) {
        auto n_nodes = x.size();
        for (std::size_t i = 0; i < n_nodes; i++) {
            bbox.xmin = std::min(bbox.xmin, x[i]);
            bbox.xmax = std::max(bbox.xmax, x[i]);
            bbox.ymin = std::min(bbox.ymin, y[i]);
            bbox.ymax = std::max(bbox.ymax, y[i]);
            bbox.zmin = std::min(bbox.zmin, z[i]);
            bbox.zmax = std::max(bbox.zmax, z[i]);
        }
        return bbox;
    }
    else
        throw std::logic_error(
            fmt::format("Size of x ({}) must equal to size of y({}) and size of z({}).",
                        x.size(),
                        y.size(),
                        z.size()));
}

/// Analyze mesh
///
/// This function analyzes the mesh and determines the dimension of the mesh
/// (1D, 2D or 3D) and the dimension of the side sets (dim - 1) and the dimension
/// of the node sets (0).
void
analyze_mesh()
{
    auto bbox = compute_mesh_bounding_box();
    double x_width = bbox.xmax - bbox.xmin;
    double y_width = bbox.ymax - bbox.ymin;
    double z_width = bbox.zmax - bbox.zmin;

    if ((std::fabs(x_width) > 0) && (std::fabs(y_width) < 1e-16) && (std::fabs(z_width) < 1e-16))
        dim = 1;
    else if ((std::fabs(x_width) > 0) && (std::fabs(y_width) > 0) && (std::fabs(z_width) < 1e-16))
        dim = 2;
    else
        dim = 3;

    side_set_dim = dim - 1;
    node_set_dim = 0;
}

/// Get physical entities with a given dimension
///
/// This function returns the physical entities with a given dimension
///
/// @param f GMSH file
/// @param dim Spatial dimension of the physical entities
/// @return Physical entities with the given dimension
const std::vector<gmshparsercpp::MshFile::MultiDEntity> &
get_entities_by_dim(const gmshparsercpp::MshFile & f, int dim)
{
    switch (dim) {
    case 1:
        return f.get_curve_entities();
    case 2:
        return f.get_surface_entities();
    case 3:
        return f.get_volume_entities();
    default:
        throw std::runtime_error(fmt::format("Unsupported spatial dimension {}.", dim));
    }
}

/// Correct the orientation of a triangle element
///
/// This function corrects the orientation of a triangle element by
/// ensuring that the cross product of the two vectors formed by the
/// edges of the triangle is positive. If the cross product is negative,
/// the order of the nodes is reversed.
///
/// @param el_nodes The nodes of the triangle element
void
correct_tri3_orientation(std::vector<int> & el_nodes)
{
    std::array<double, 2> pt1 = { x[el_nodes[0] - 1], y[el_nodes[0] - 1] };
    std::array<double, 2> pt2 = { x[el_nodes[1] - 1], y[el_nodes[1] - 1] };
    std::array<double, 2> pt3 = { x[el_nodes[2] - 1], y[el_nodes[2] - 1] };

    std::array<double, 3> vec1 = { pt2[0] - pt1[0], pt2[1] - pt1[1] };
    std::array<double, 3> vec2 = { pt3[0] - pt1[0], pt3[1] - pt1[1] };

    double n = vec1[0] * vec2[1] - vec1[1] * vec2[0];
    if (n < 0) {
        // this should never happen, since gmsh claims to produce CCW oriented triangles, but it
        // does *sigh* so, reorient the triangle...
        std::swap(el_nodes[1], el_nodes[2]);
    }
}

/// Build coordinates for the exodus file
///
/// This function builds the coordinates for the exodus file. It also builds the node map, which
/// maps the GMSH node id to internal node id.
///
/// @param nodes The nodes from the GMSH file
void
build_coordinates(const std::vector<gmshparsercpp::MshFile::Node> & nodes)
{
    for (const auto & nd : nodes) {
        if (!nd.tags.empty()) {
            for (int j = 0; j < nd.tags.size(); j++) {
                auto local_id = (int) node_map.size();
                const auto & id = nd.tags[j];
                node_map[id] = local_id;

                const auto & c = nd.coordinates[j];
                x.push_back(c.x);
                y.push_back(c.y);
                z.push_back(c.z);
            }
        }
    }
}

/// Build element block by dimension
///
/// This function builds the element blocks by dimension. It fills the `el_blk_dim` vector with the
/// element blocks based on their dimension.
///
/// @param el_blks The element blocks from the GMSH file
void
build_element_block_dim(const std::vector<gmshparsercpp::MshFile::ElementBlock> & el_blks)
{
    el_blk_dim.resize(4);
    for (const auto & eb : el_blks)
        el_blk_dim[eb.dimension].push_back(&eb);
}

void
build_element_blocks(const std::vector<const gmshparsercpp::MshFile::ElementBlock *> & el_blks,
                     const std::vector<gmshparsercpp::MshFile::MultiDEntity> & entities)
{
    std::map<int, const gmshparsercpp::MshFile::MultiDEntity *> ents_by_id;
    for (const auto & ent : entities)
        ents_by_id[ent.tag] = &ent;

    std::map<int, std::string> exo_names;
    std::map<int, std::vector<int>> exo_connects;
    std::map<int, gmshparsercpp::ElementType> eltype;
    std::map<int, int> num_elems;

    unsigned int eid = 0;
    for (const auto & eb : el_blks) {
        const auto * ent = ents_by_id[eb->tag];
        int id;
        if (ent == nullptr)
            id = eb->tag;
        else if (!ent->physical_tags.empty()) {
            // the sign on physical tag ID refers to orientation which we don't need
            id = std::abs(ent->physical_tags[0]);
            auto it = phys_ent_by_tag.find(id);
            if (it != phys_ent_by_tag.end())
                exo_names[id] = it->second->name;
        }
        else
            id = eb->tag;

        if (eltype.find(id) == eltype.end())
            eltype[id] = eb->element_type;
        else if (eltype[id] != eb->element_type) {
            fmt::print("{} {} | {}\n", eltype[id], eb->element_type, gmshparsercpp::NONE);
            throw std::runtime_error(
                fmt::format("Unable to combine different element types in a single block"));
        }

        std::vector<int> & connect = exo_connects[id];
        for (const auto & elem : eb->elements) {
            std::vector<int> el_nodes;
            auto n_node_tags = elem.node_tags.size();
            for (std::size_t i = 0; i < n_node_tags; i++) {
                const auto & nid = elem.node_tags[node_order[eb->element_type][i]];
                el_nodes.push_back(nid);
            }
            if (eb->element_type == gmshparsercpp::TRI3)
                correct_tri3_orientation(el_nodes);
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
                // See [1] for how sides are numbered on different elements (fig 4.15, pp. 28)
                if (eb->element_type == gmshparsercpp::TET4) {
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
                else if (eb->element_type == gmshparsercpp::HEX8) {
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

        num_elems[id] += (int) eb->elements.size();
    }

    for (auto const & [id, connect] : exo_connects) {
        exodusIIcpp::ElementBlock eb;

        eb.set_name(exo_names[id]);
        eb.set_id(id);
        auto elem_type = eltype[id];
        eb.set_connectivity(exo_elem_type.at(elem_type),
                            num_elems[id],
                            nodes_per_elem[elem_type],
                            connect);
        element_blocks.push_back(eb);
    }
}

/// Build (0-D) side sets from the point entities
///
/// @param el_blks Element blocks
/// @param entities Point entities
void
build_side_sets(const std::vector<const gmshparsercpp::MshFile::ElementBlock *> & el_blks,
                const std::vector<gmshparsercpp::MshFile::PointEntity> & entities)
{
    std::map<int, const gmshparsercpp::MshFile::PointEntity *> ents_by_id;
    for (const auto & ent : entities) {
        if (!ent.physical_tags.empty()) {
            for (const auto & tag : ent.physical_tags) {
                exodusIIcpp::SideSet ss;
                // the sign on physical tag ID refers to orientation which we don't need
                auto id = std::abs(tag);
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
        for (const auto & tag : ent->physical_tags) {
            // the sign on physical tag ID refers to orientation which we don't need
            auto id = std::abs(tag);
            auto & ss = side_sets[id];
            for (const auto & elem : eb->elements) {
                if (eb->element_type == gmshparsercpp::POINT) {
                    std::vector<int> side_key = { elem.node_tags[0] };
                    const auto el_side_pair = elem_sides[side_key];
                    ss.add(el_side_pair.first, el_side_pair.second);
                }
            }
        }
    }
}

/// Build (1-D, 2-D) side sets from the physical entities
///
/// @param el_blks Element blocks
/// @param entities Physical entities
void
build_side_sets(const std::vector<const gmshparsercpp::MshFile::ElementBlock *> & el_blks,
                const std::vector<gmshparsercpp::MshFile::MultiDEntity> & entities)
{
    std::map<int, const gmshparsercpp::MshFile::MultiDEntity *> ents_by_id;
    for (const auto & ent : entities) {
        if (!ent.physical_tags.empty()) {
            for (const auto & tag : ent.physical_tags) {
                exodusIIcpp::SideSet ss;
                // the sign on physical tag ID refers to orientation which we don't need
                auto id = std::abs(tag);
                auto it = phys_ent_by_tag.find(id);
                if (it != phys_ent_by_tag.end())
                    ss.set_name(it->second->name);
                side_sets[id] = ss;
            }
        }
        ents_by_id[ent.tag] = &ent;
    }

    for (const auto & eb : el_blks) {
        int id;
        if (ents_by_id.count(eb->tag) > 0) {
            const auto * ent = ents_by_id[eb->tag];
            if (ent->physical_tags.size() > 0)
                id = std::abs(ent->physical_tags[0]);
            else
                id = eb->tag;
        }
        else
            id = eb->tag;

        auto & ss = side_sets[id];
        ss.set_id(id);
        for (const auto & elem : eb->elements) {
            std::vector<int> side_key;
            switch (eb->element_type) {
            case gmshparsercpp::LINE2:
                side_key = build_side_key_edge2(elem.node_tags[0], elem.node_tags[1]);
                break;
            case gmshparsercpp::TRI3:
                side_key =
                    build_side_key_tri3(elem.node_tags[0], elem.node_tags[1], elem.node_tags[2]);
                break;
            case gmshparsercpp::QUAD4:
                side_key = build_side_key_quad4(elem.node_tags[0],
                                                elem.node_tags[1],
                                                elem.node_tags[2],
                                                elem.node_tags[3]);
                break;
            default:
                break;
            }
            const auto el_side_pair = elem_sides[side_key];
            ss.add(el_side_pair.first, el_side_pair.second);
        }
    }
}

///
void
build_node_sets(const std::vector<const gmshparsercpp::MshFile::ElementBlock *> & el_blks,
                const std::vector<gmshparsercpp::MshFile::PointEntity> & entities)
{
    std::map<int, const gmshparsercpp::MshFile::PointEntity *> ents_by_id;
    for (const auto & ent : entities) {
        if (!ent.physical_tags.empty()) {
            for (const auto & tag : ent.physical_tags) {
                exodusIIcpp::NodeSet ns;
                // the sign on physical tag ID refers to orientation which we don't need
                auto id = std::abs(tag);
                ns.set_id(id);
                auto it = phys_ent_by_tag.find(id);
                if (it != phys_ent_by_tag.end())
                    ns.set_name(it->second->name);
                node_sets[id] = ns;
            }
        }
        ents_by_id[ent.tag] = &ent;
    }

    for (const auto & eb : el_blks) {
        const auto & ent = ents_by_id[eb->tag];
        for (const auto & tag : ent->physical_tags) {
            std::vector<int> nodes;
            for (const auto & elem : eb->elements) {
                if (eb->element_type == gmshparsercpp::POINT) {
                    auto node_key = elem.node_tags[0];
                    auto node_id = node_map[node_key];
                    nodes.push_back(node_id + 1);
                }
            }

            // the sign on physical tag ID refers to orientation which we don't need
            auto id = std::abs(tag);
            node_sets[id].set_nodes(nodes);
        }
    }
}

/// Read GMSH file
///
/// @param file_name GMSH file name
void
read_gmsh_file(const std::string & file_name)
{
    gmshparsercpp::MshFile f(file_name);
    f.parse();

    const auto & phys_entities = f.get_physical_names();
    read_physical_entities(phys_entities);

    const auto & nodes = f.get_nodes();
    build_coordinates(nodes);
    const auto & el_blks = f.get_element_blocks();
    analyze_mesh();
    build_element_block_dim(el_blks);

    const auto & ents = get_entities_by_dim(f, dim);
    build_element_blocks(el_blk_dim[dim], ents);

    if (side_set_dim == 0) {
        const auto & sideset_ents = f.get_point_entities();
        build_side_sets(el_blk_dim[side_set_dim], sideset_ents);
    }
    else {
        const auto & sideset_ents = get_entities_by_dim(f, side_set_dim);
        build_side_sets(el_blk_dim[side_set_dim], sideset_ents);
    }

    const auto & nodeset_ents = f.get_point_entities();
    build_node_sets(el_blk_dim[node_set_dim], nodeset_ents);
}

/// Write coordinates to ExodusII file
///
/// @param f ExodusII file
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

/// Write element blocks to ExodusII file
///
/// @param f ExodusII file
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

/// Write side sets to ExodusII file
///
/// @param f ExodusII file
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

/// Write node sets to ExodusII file
///
/// @param f ExodusII file
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

/// Write ExodusII file
///
/// @param file_name ExodusII file name
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

/// Convert GMSH file to ExodusII file
///
/// @param input_file_name GMSH file name
/// @param output_file_name ExodusII file name
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
        if (result.count("version")) {
            fmt::print(stdout, "gmsh2exo version {}\n", GMSH2EXO_VERSION);
        }
        else if (result.count("gmsh_file") && result.count("exo2_file")) {
            convert(result["gmsh_file"].as<std::string>(), result["exo2_file"].as<std::string>());
        }
        else {
            fmt::print(stdout, options.help());
        }
    }
    catch (const cxxopts::exceptions::exception & e) {
        fmt::print(stderr, "Error: {}\n", e.what());
        fmt::print(stdout, options.help());
        return 1;
    }
    catch (std::exception & e) {
        fmt::print(stderr, "Error: {}\n", e.what());
        return 1;
    }

    return 0;
}

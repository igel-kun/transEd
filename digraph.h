/***************************************************
 * digraph.h
  * by M.Weller
 **************************************************/
 
#ifndef digraph_h
#define digraph_h

#include "stdlib.h"
#include "string.h"
#include <set> 
#include <map>
#include <string>
#include <vector>

#define DEBUG 0
#define dbgcout if(DEBUG) cout

#define MY_INFINITY (uint)-2
#define MRK_NONE 0
#define MRK_FORBIDDEN 1
#define MRK_PERMANENT 2

using namespace std;

typedef unsigned int uint;

typedef string t_vertex;
typedef pair<t_vertex, t_vertex> t_arc;
typedef pair<t_arc, uint> t_labeled_arc; // 0 = no label, 1 = forbidden, 2 = permanent
typedef set<t_vertex> t_adjlist;
typedef map<t_vertex, t_adjlist> t_adjtable;
typedef pair<t_arc, t_arc> t_3path;
typedef t_arc t_diamond; // only the first and the last vertex of a diamond 
						// are saved, the belt can be calced in O(n)

typedef t_3path* p_3path;
typedef t_diamond* p_diamond;

#include <iostream>
#include <fstream>
#include <algorithm>

inline set<t_arc> get_irreflexive(const set<t_arc>& s){
	set<t_arc> result;
	for(set<t_arc>::const_iterator i = s.begin(); i != s.end(); i++)
		if(i->first != i->second) result.insert(*i);
	return result;
}
inline set<t_vertex> set_unite(const set<t_vertex>& s1, const set<t_vertex>& s2){
	set<t_vertex> result(s1);
//	set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), result.begin());
//	for(set<t_vertex>::const_iterator i = s2.begin(); i != s2.end(); i++) result.insert(*i);
	result.insert(s2.begin(), s2.end());
	return result;
}
inline set<t_vertex> set_intersect(const set<t_vertex>& s1, const set<t_vertex>& s2){
	set<t_vertex> result;
//	for(set<t_vertex>::const_iterator i = s1.begin(); i != s1.end(); i++)
//		if(s2.find(*i) != s2.end()) result.insert(*i);
	set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(result,result.begin()));
	return result;
}
inline set<t_vertex> set_substract(const set<t_vertex>& s1, const set<t_vertex>& s2){
//	set<t_vertex> result(s1);
//	for(set<t_vertex>::const_iterator i = s2.begin(); i != s2.end(); i++) result.erase(*i);
	set<t_vertex> result;
	set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(result,result.begin()));
	return result;
}
inline set<t_vertex> set_symm_diff(const set<t_vertex>& s1, const set<t_vertex>& s2){
//	set<t_vertex> result(s1);
//	for(set<t_vertex>::const_iterator i = s2.begin(); i != s2.end(); i++) 
//		if(s2.find(*i) != s2.end()) result.erase(*i); else result.insert(*i);
	set<t_vertex> result;
	set_symmetric_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(result,result.begin()));
	return result;
}
inline set<t_arc> set_multiply(const set<t_vertex>& s1, const set<t_vertex>& s2){
	set<t_arc> result;
	for(set<t_vertex>::const_iterator i = s1.begin(); i != s1.end(); i++) 
		for(set<t_vertex>::const_iterator j = s2.begin(); j != s2.end(); j++)
			result.insert(t_arc(*i,*j));
	return get_irreflexive(result);
}
inline set<t_arc> set_substract(const set<t_arc>& s1, const set<t_arc>& s2){
//	set<t_arc> result(s1);
//	for(set<t_arc>::const_iterator i = s2.begin(); i != s2.end(); i++) result.erase(*i);
	set<t_arc> result;
	set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(result,result.begin()));
	return result;
}

/*
set<t_vertex> set_unite(const set<t_vertex>& s1, const set<t_vertex>& s2);
set<t_vertex> set_intersect(const set<t_vertex>& s1, const set<t_vertex>& s2);
set<t_vertex> set_substract(const set<t_vertex>& s1, const set<t_vertex>& s2);
set<t_vertex> set_symm_diff(const set<t_vertex>& s1, const set<t_vertex>& s2);
set<t_arc> set_multiply(const set<t_vertex>& s1, const set<t_vertex>& s2);
set<t_arc> set_substract(const set<t_arc>& s1, const set<t_arc>& s2);
set<t_arc> get_irreflexive(const set<t_arc>& s);
*/
inline t_arc reverse_arc(const t_arc a){
	return t_arc(a.second, a.first);
}

set<t_vertex> get_maximal_matching(const set<t_vertex>& vertices, const set<t_arc>& arcs);

ostream& operator<<(ostream& os, const set<t_vertex>& s);
ostream& operator<<(ostream& os, const set<t_arc>& a);

ostream& operator<<(ostream& os, const t_arc& a);

class t_dir_graph {
private:
	t_adjtable successors;
	t_adjtable predecessors;
	set<t_arc> forbidden_arcs;
	set<t_arc> permanent_arcs;

	set<t_3path> p3_set;
	map<t_diamond, uint> diamond_set; // map diamonds to belt sizes
public:
// ============ construction / destruction =============
	// constructor
	t_dir_graph();
	// copy constructor
	t_dir_graph(const t_dir_graph& d);
	// destructor
	~t_dir_graph();
// ============== operators ==================
	t_dir_graph operator=(const t_dir_graph& d);
// ============== information providing ==============
	// return whether the digraph is empty [contains no arcs]
	bool has_no_arcs() const;
	// return whether the digraph is empty [contains no vertices]
	bool is_empty() const;
	// return whether the digraph is completely connected
	bool is_complete() const;
	// return whether the given vertex is in the digraph
	bool contains_vertex(const t_vertex& v) const;
	// return whether the given arc is in the digraph
	bool contains_arc(const t_arc& a) const;
	// return whether the given arc is forbidden
	bool is_forbidden(const t_arc& a) const;
	// return whether the given arc is permanent
	bool is_permanent(const t_arc& a) const;
	// return whether v is a sink of the digraph
	bool is_sink(const t_vertex& v) const;
	// return whether v is a source of the digraph
	bool is_source(const t_vertex& v) const;
	// return all sinks
	set<t_vertex> get_sinks() const;
	// return all sources
	set<t_vertex> get_sources() const;
	// get all vertices
	set<t_vertex> get_vertices() const;
	// get all arcs
	set<t_arc> get_arcs() const;
	// get all successor-vertices of a given vertex
	set<t_vertex> get_successors(const t_vertex& v) const;
	// get all predecessor-vertices of a given vertex
	set<t_vertex> get_predecessors(const t_vertex& v) const;
	// get all neighbors of v, that is successors and predecessors
	set<t_vertex> get_neighbors(const t_vertex& v) const;
	// get the number of vertices in the graph
	uint get_vertex_number() const;
	// get the number of arcs in the graph
	uint get_arc_number() const;
	// get the subgraph induced by the given vertex-set
	t_dir_graph* get_induced_subgraph(const set<t_vertex>& vertices) const;
	// calculate the reachability of v in D
	// if undirected is true, all vertices can reach their predecessors
	set<t_vertex> get_reachable(const t_vertex& v, bool undirected = false) const;
// ========== P3 information providing ==============
	// return all p3 in p3_set with (v,*,*)
	set<t_3path> get_p3_by_first(const t_vertex& v) const;
	// return all p3 in p3_set with (*,v,*)
	set<t_3path> get_p3_by_middle(const t_vertex& v) const;
	// return all p3 in p3_set with (*,*,v)
	set<t_3path> get_p3_by_third(const t_vertex& v) const;
	// get the set of P3 that involve a given vertex v
	set<t_3path> get_p3_involving(const t_vertex& v) const;

	// get some p3 in the digraph
	p_3path get_a_p3() const;
	// get some diamond in the digraph
	p_diamond get_a_diamond() const;

	// get all p3 in the digraph
	set<t_3path> get_p3() const;
	// get all diamonds in the digraph
	map<t_diamond,uint> get_diamonds() const;
	// get all vertices that are the first(second,third) vertex of some P3
	set<t_vertex> get_p3_first() const;
	set<t_vertex> get_p3_second() const;
	set<t_vertex> get_p3_third() const;
// ========== graph modifications ===================
	// increase the belt size of the given diamond
	void increase_diamond(const t_diamond& diam);
	// decrease the belt size of the given diamond
	void decrease_diamond(const t_diamond& diam);
	// insert a vertex
	bool insert_vertex(const t_vertex& v);
	// delete a vertex
	bool delete_vertex(const t_vertex& v);
	// delete multiple vertices
	bool delete_vertices(const set<t_vertex>& V);
	// insert an arc
	bool insert_arc(const t_arc& a, const bool permanent = false);
	// delete an arc
	bool delete_arc(const t_arc& a, const bool forbidden = false);
	// mark an arc forbidden/permanent
	bool mark_arc(const t_arc& a, const uint mark);
	// delete all vertices that are not in the given vertex set from the graph
	void induced_subgraph(const set<t_vertex>& vertices);
	// split a weakly connected component from the digraph
	t_dir_graph* split_component();
// ============= IO ====================
	// read a graph from the given file, return success
	bool read_from_stream(istream& input);
	// read a graph from the given file, return success
	bool read_from_file(const char* filename);
	// print the digraph to a stream
	void print(ostream& output) const;
	// streaming IO
	friend ostream& operator<<(ostream& os, const t_dir_graph& d);
	friend istream& operator>>(istream& is, t_dir_graph& d);
private:
	bool read_format1(istream& input, const string first_line);
	bool read_format2(istream& input, const string first_line);
};


t_dir_graph* generate_random_digraph(const uint vertices, const double arc_probability = 0.5);


#endif

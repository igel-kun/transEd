/***************************************************
 * digraph.cpp
  * by M.Weller
 **************************************************/
#ifndef digraph_cpp
#define digraph_cpp

#include "digraph.h"




// return a maximal matching of the digraph (obtained greedily)
set<t_vertex> get_maximal_matching(const set<t_vertex>& vertices, const set<t_arc>& arcs){
	set<t_vertex> result;
	set<t_arc> arcs_left = arcs;
	set<t_arc> to_remove;
	t_vertex u,v;

	while(!arcs_left.empty()){
		u = arcs_left.begin()->first;
		v = arcs_left.begin()->second;

		result.insert(u);
		result.insert(v);

		to_remove.clear();
		for(set<t_arc>::iterator a = arcs_left.begin(); a != arcs_left.end(); a++)
			if(((a->first == u)) || (a->first == v) || (a->second == u) || (a->second == v)) 
				to_remove.insert(*a);
		arcs_left = set_substract(arcs_left, to_remove);
	}
	return result;
}

ostream& operator<<(ostream& os, const set<t_vertex>& s){
	if(s.size()) {
		os << "(";
		for(set<t_vertex>::const_iterator i = s.begin(); i != s.end(); i++)
			os << *i <<",";
		os << "\b)";
	} else os << "()";
	return os;
}

ostream& operator<<(ostream& os, const set<t_arc>& a){
	if(a.size()){
		os << "(";
		for(set<t_arc>::const_iterator i = a.begin(); i != a.end(); i++)
			os << *i << ",";
		os << "\b)";
	} else os << "()";
	return os;
}	
// CLASS: t_dir_graph
// ============ construction / destruction =============
// constructor
t_dir_graph::t_dir_graph(){}
// copy constructor
t_dir_graph::t_dir_graph(const t_dir_graph& d)
	:successors(d.successors), predecessors(d.predecessors), forbidden_arcs(d.forbidden_arcs), permanent_arcs(d.permanent_arcs), p3_set(d.p3_set), diamond_set(d.diamond_set){}
// destructor
t_dir_graph::~t_dir_graph(){}
t_dir_graph t_dir_graph::operator=(const t_dir_graph& d){
	successors = d.successors;
	predecessors = d.predecessors;
	forbidden_arcs = d.forbidden_arcs;
	permanent_arcs = d.permanent_arcs;
	p3_set = d.p3_set;
	diamond_set = d.diamond_set;
	return *this;
}
// ============== information providing ==============
// return whether the digraph is empty [contains no arcs]
bool t_dir_graph::has_no_arcs() const{
	for(t_adjtable::const_iterator i = successors.begin(); i != successors.end(); i++)
		if(!i->second.empty()) return false;
	return true;
}

// return whether the digraph is empty [contains no vertices]
bool t_dir_graph::is_empty() const{
	return (successors.size() == 0);
}

// return whether the digraph is completely connected
bool t_dir_graph::is_complete() const{
	for(t_adjtable::const_iterator i = successors.begin(); i != successors.end(); i++)
		if(i->second.size() != successors.size() - 1) return false;
	return true;
}
// return whether the given vertex is in the digraph
bool t_dir_graph::contains_vertex(const t_vertex& v) const{
	return (successors.find(v) != successors.end());
}
// return whether the given arc is in the digraph
bool t_dir_graph::contains_arc(const t_arc& a) const {
	t_adjtable::const_iterator i = successors.find(a.first);
	if(i != successors.end()) {
		return i->second.find(a.second) != i->second.end();
	} else return false;
}
// return whether the given arc is forbidden
bool t_dir_graph::is_forbidden(const t_arc& a) const{
	return (forbidden_arcs.find(a) != forbidden_arcs.end());
}
// return whether the given arc is permanent
bool t_dir_graph::is_permanent(const t_arc& a) const{
	return (permanent_arcs.find(a) != permanent_arcs.end());
}

// return whether v is a sink of the digraph
bool t_dir_graph::is_sink(const t_vertex& v) const{
	t_adjtable::const_iterator i = successors.find(v);
	if(i != successors.end()) return i->second.empty(); else return false;
}
// return whether v is a source of the digraph
bool t_dir_graph::is_source(const t_vertex& v) const{
	t_adjtable::const_iterator i = predecessors.find(v);
	if(i != predecessors.end()) return i->second.empty(); else return false;
}
// return all sinks
set<t_vertex> t_dir_graph::get_sinks() const{
	set<t_vertex> result;
	for(t_adjtable::const_iterator i = successors.begin(); i != successors.end(); i++)
		if(i->second.empty()) result.insert(i->first);
	return result;
}
// return all sources
set<t_vertex> t_dir_graph::get_sources() const{
	set<t_vertex> result;
	for(t_adjtable::const_iterator i = predecessors.begin(); i != predecessors.end(); i++)
		if(i->second.empty()) result.insert(i->first);
	return result;
}
// get all vertices
set<t_vertex> t_dir_graph::get_vertices() const{
	set<t_vertex> result;
	for(t_adjtable::const_iterator i = successors.begin(); i != successors.end(); i++)
		result.insert(i->first);
	return result;
}
// get all arcs
set<t_arc> t_dir_graph::get_arcs() const{
	set<t_arc> result;
	for(t_adjtable::const_iterator i = successors.begin(); i != successors.end(); i++)
		for(t_adjlist::const_iterator j = i->second.begin(); j != i->second.end(); j++)
			result.insert(t_arc(i->first,*j));
	return result;
}
// get all adjacent vertices of a given vertex
set<t_vertex> t_dir_graph::get_successors(const t_vertex& v) const{
	t_adjtable::const_iterator i = successors.find(v);
	if(i != successors.end())
		return i->second;
	else return set<t_vertex>();
}
// get all predecessor-vertices of a given vertex
set<t_vertex> t_dir_graph::get_predecessors(const t_vertex& v) const{
	t_adjtable::const_iterator i = predecessors.find(v);
	if(i != predecessors.end())
		return i->second;
	else return set<t_vertex>();
}
// get all neighbors of v, that is successors and predecessors
set<t_vertex> t_dir_graph::get_neighbors(const t_vertex& v) const{
	return set_unite(get_predecessors(v), get_successors(v));
}
// get the number of vertices in the graph
uint t_dir_graph::get_vertex_number() const{
	return successors.size();
}
// get the number of arcs in the graph
uint t_dir_graph::get_arc_number() const{
	uint result = 0;
	for(t_adjtable::const_iterator i = successors.begin(); i != successors.end(); i++)
		result += i->second.size();
	return result;
}
// get the subgraph induced by the given vertex-set
t_dir_graph* t_dir_graph::get_induced_subgraph(const set<t_vertex>& vertices) const{
	t_dir_graph* result = new t_dir_graph(*this);
	result->induced_subgraph(vertices);
	return result;
}
// calculate the reachability of v in D
// if undirected is true, all vertices can reach their predecessors
set<t_vertex> t_dir_graph::get_reachable(const t_vertex& v, bool undirected) const{
	if(contains_vertex(v)){
		set<t_vertex> to_check, checked, neighbors;
		t_vertex u;

		to_check.insert(v);
		while(!to_check.empty()){
			u = *(to_check.begin());
			to_check.erase(to_check.begin());
			
			neighbors = get_successors(u);
			if(undirected) neighbors = set_unite(neighbors, get_predecessors(u));
			
			neighbors = set_substract(neighbors, checked);
			to_check = set_unite(to_check, neighbors);

			checked.insert(u);
		}
		return checked;
	} else return set<t_vertex>();
}

// return all p3 in p3_set with (v,*,*)
set<t_3path> t_dir_graph::get_p3_by_first(const t_vertex& v) const{
	set<t_vertex> succ, succsucc;
	set<t_3path> result;
	succ = get_successors(v);
	for(set<t_vertex>::const_iterator u = succ.begin(); u!= succ.end(); u++){
		succsucc = get_successors(*u);
		for(set<t_vertex>::const_iterator w = succsucc.begin(); w!= succsucc.end(); w++)
			if(p3_set.find(t_3path(t_arc(v,*u),t_arc(*u,*w))) != p3_set.end())
				result.insert(t_3path(t_arc(v,*u),t_arc(*u,*w)));
	}
	return result;
}

// return all p3 in p3_set with (*,v,*)
set<t_3path> t_dir_graph::get_p3_by_middle(const t_vertex& v) const{
	set<t_vertex> succ, pred;
	set<t_3path> result;
	succ = get_successors(v);
	pred = get_predecessors(v);
	for(set<t_vertex>::const_iterator u = pred.begin(); u!= pred.end(); u++)
		for(set<t_vertex>::const_iterator w = succ.begin(); w!= succ.end(); w++)
			if(p3_set.find(t_3path(t_arc(*u,v),t_arc(v,*w))) != p3_set.end())
				result.insert(t_3path(t_arc(*u,v),t_arc(v,*w)));
	return result;
}

// return all p3 in p3_set with (*,*,v)
set<t_3path> t_dir_graph::get_p3_by_third(const t_vertex& v) const{
	set<t_vertex> pred, predpred;
	set<t_3path> result;
	pred = get_predecessors(v);
	for(set<t_vertex>::const_iterator u = pred.begin(); u!= pred.end(); u++){
		predpred = get_predecessors(*u);
		for(set<t_vertex>::const_iterator w = predpred.begin(); w!= predpred.end(); w++)
			if(p3_set.find(t_3path(t_arc(*w,*u),t_arc(*u,v))) != p3_set.end())
				result.insert(t_3path(t_arc(*w,*u),t_arc(*u,v)));
	}
	return result;
}
// get the set of P3 that involve a given vertex v
// O(n^2 log n)
set<t_3path> t_dir_graph::get_p3_involving(const t_vertex& v) const{
	set<t_3path> result = get_p3_by_first(v);
	set<t_3path> tmp = get_p3_by_middle(v);
	result.insert(tmp.begin(), tmp.end());
	tmp = get_p3_by_third(v);
	result.insert(tmp.begin(), tmp.end());
	return result;
}

// get some p3 in the digraph
p_3path t_dir_graph::get_a_p3() const{
	if(p3_set.empty()) return NULL;
		else return new t_3path(*p3_set.begin());
}
// get some diamond in the digraph
p_diamond t_dir_graph::get_a_diamond() const{
	if(diamond_set.empty()) return NULL;
		else return new t_diamond(diamond_set.begin()->first);
}

// get all p3 in the digraph
set<t_3path> t_dir_graph::get_p3() const{
	return p3_set;
}
// get all diamonds in the digraph
map<t_diamond,uint> t_dir_graph::get_diamonds() const{
	return diamond_set;
}

// get all vertices that are the first(second,third) vertex of some P3
set<t_vertex> t_dir_graph::get_p3_first() const{
	set<t_vertex> result;
	for(set<t_3path>::const_iterator u = p3_set.begin(); u!= p3_set.end(); u++)
		result.insert(u->first.first);
	return result;
}
set<t_vertex> t_dir_graph::get_p3_second() const{
	set<t_vertex> result;
	for(set<t_3path>::const_iterator u = p3_set.begin(); u!= p3_set.end(); u++)
		result.insert(u->first.second);
	return result;
}
set<t_vertex> t_dir_graph::get_p3_third() const{
	set<t_vertex> result;
	for(set<t_3path>::const_iterator u = p3_set.begin(); u!= p3_set.end(); u++)
		result.insert(u->second.second);
	return result;
}


// ========== graph modifications ===================

// increase the belt size of the given diamond
// O(log n)
void t_dir_graph::increase_diamond(const t_diamond& diam){
	map<t_diamond, uint>::iterator d_iter = diamond_set.find(diam);
//	cout << "increasing diamond-belt of "<<diam<<"\n";
	if(d_iter != diamond_set.end()){
		d_iter->second++;
	} else {// if its not there, maybe it should be inserted...
		// this intersection can be calced in O(1), otherwise, diam would have been in diamond_set !
		set<t_vertex> belt = set_intersect(get_successors(diam.first), get_predecessors(diam.second));
//		cout << "oops, not yet there, "<<diam<<" is "<<(belt.size()>1?"":"not ")<<"inserted (size "<<belt.size()<<")\n";
		if(belt.size() > 1)
			diamond_set.insert(pair<t_diamond,uint>(diam, belt.size()));
	}
}
// decrease the belt size of the given diamond
// O(log n)
void t_dir_graph::decrease_diamond(const t_diamond& diam){
	map<t_diamond, uint>::iterator d_iter = diamond_set.find(diam);
	if(d_iter != diamond_set.end()){
//		cout << "decreasing diamond-belt of "<<diam<<" (current size: "<<d_iter->second <<")\n";
		if(d_iter->second > 2) d_iter->second--;
			else diamond_set.erase(d_iter);
	}
}


// insert a vertex
// O(log n)
bool t_dir_graph::insert_vertex(const t_vertex& v){
	pair<t_adjtable::iterator, bool> result;
	result = successors.insert(pair<t_vertex,t_adjlist>(v,t_adjlist()));
	if(result.second){
		if(predecessors.insert(pair<t_vertex,t_adjlist>(v,t_adjlist())).second)
			return true;
		else {// if insertion failed, remove the vertex from the succ-list as well to avoid inconsistency
			successors.erase(result.first);
			return false;
		}
	} else return false;
}
// delete a vertex
// O(n^2)
bool t_dir_graph::delete_vertex(const t_vertex& v){
	// 1. update the P3-set & the diamond-set
	set<t_3path> bad_p3 = get_p3_involving(v);

	for(set<t_3path>::const_iterator i = bad_p3.begin(); i != bad_p3.end(); i++)
		decrease_diamond(t_diamond(i->first.first, i->second.second));

	set<t_3path> temp_p3;
	set_difference(p3_set.begin(), p3_set.end(), bad_p3.begin(), bad_p3.end(), inserter(temp_p3,temp_p3.begin()));
	p3_set = temp_p3;
	// 2. delete all arcs that involve v
	for(t_adjtable::iterator i = successors.begin(); i != successors.end(); i++)
		i->second.erase(v);
	for(t_adjtable::iterator i = predecessors.begin(); i != predecessors.end(); i++)
		i->second.erase(v);
	// 3. delete v itself
	return (successors.erase(v)>0 && (predecessors.erase(v)>0));
}
// delete multiple vertices
bool t_dir_graph::delete_vertices(const set<t_vertex>& V){
	// 1. update the P3-set & the diamond-set
	set<t_3path> bad_p3, new_bad_p3;
	for(set<t_vertex>::const_iterator i = V.begin(); i != V.end(); i++){
		new_bad_p3 = get_p3_involving(*i);
		bad_p3.insert(new_bad_p3.begin(), new_bad_p3.end());
	}
	for(set<t_3path>::const_iterator i = bad_p3.begin(); i != bad_p3.end(); i++)
		decrease_diamond(t_diamond(i->first.first, i->second.second));

	set<t_3path> temp_p3;
	set_difference(p3_set.begin(), p3_set.end(), bad_p3.begin(), bad_p3.end(), inserter(temp_p3,temp_p3.begin()));
	p3_set = temp_p3;
	// 2. delete all arcs that involve V
	for(t_adjtable::iterator i = successors.begin(); i != successors.end(); i++)
		i->second = set_substract(i->second, V);
	for(t_adjtable::iterator i = predecessors.begin(); i != predecessors.end(); i++)
		i->second = set_substract(i->second, V);
	// 3. delete V itself
	bool result = true;
	for(set<t_vertex>::iterator i = V.begin(); i != V.end(); i++){
		result &= (successors.erase(*i)>0);
		result &= (predecessors.erase(*i)>0);
	}
	return result;
}
// insert an arc
// O(n log n)
bool t_dir_graph::insert_arc(const t_arc& a, const bool permanent){
	if(successors.find(a.second) == successors.end()) return false; // dont insert edges that are not in VxV
	if(predecessors.find(a.first) == predecessors.end()) return false; // dont insert edges that are not in VxV
	if(forbidden_arcs.find(a) != forbidden_arcs.end()) return false; // dont insert if the arc is forbidden
	t_adjtable::iterator i = successors.find(a.first);
	t_adjtable::iterator j = predecessors.find(a.second);
	if((i != successors.end()) && (j != predecessors.end()))
		if(i->second.insert(a.second).second)
			if(j->second.insert(a.first).second){
				if(permanent) mark_arc(a, MRK_PERMANENT);

				// update p3_set and diamond_set
				// first, remove all p3 destroyed by a
				set<t_vertex> belt = set_intersect(i->second,j->second);
				for(set<t_vertex>::const_iterator u = belt.begin(); u != belt.end(); u++)
					p3_set.erase(t_3path(t_arc(a.first, *u),t_arc(*u, a.second)));
				diamond_set.erase(t_diamond(a.first, a.second));

				// second, add new P3 and diamonds
				// u -> a.first -> a.second
				set<t_vertex> pred = set_substract(get_predecessors(a.first), j->second);
				pred.erase(a.second);
				for(set<t_vertex>::const_iterator u = pred.begin(); u != pred.end(); u++){
					increase_diamond(t_diamond(*u, a.second));
					p3_set.insert(t_3path(t_arc(*u,a.first),a));
				}

				// a.first -> a.second -> u
				set<t_vertex> succ = set_substract(get_successors(a.second), i->second);
				succ.erase(a.first);
				for(set<t_vertex>::const_iterator u = succ.begin(); u != succ.end(); u++){
					increase_diamond(t_diamond(a.first, *u));
					p3_set.insert(t_3path(a,t_arc(a.second, *u)));
				}
				return true;
			} else {// avoid inconsistency
				i->second.erase(a.second);
				return false;
			}
		else return false;
	else return false;
}
// delete an arc
// O(n log n)
bool t_dir_graph::delete_arc(const t_arc& a, const bool forbidden){
	if(permanent_arcs.find(a) != permanent_arcs.end()) return false; // dont delete if the arc is permanent
	t_adjtable::iterator i = successors.find(a.first);
	t_adjtable::iterator j = predecessors.find(a.second);
	if((i != successors.end()) && (j != predecessors.end()))
		if(i->second.erase(a.second))
			if(j->second.erase(a.first)){
				if(forbidden) mark_arc(a, MRK_FORBIDDEN);
				// update p3_set and diamond_set

				// first, remove all p3 destroyed by a
				// u -> a.first -> a.second
				set<t_vertex> pred = set_substract(get_predecessors(a.first), j->second);
				pred.erase(a.second);
				for(set<t_vertex>::const_iterator u = pred.begin(); u != pred.end(); u++){
					decrease_diamond(t_diamond(*u,a.second));
					p3_set.erase(t_3path(t_arc(*u,a.first),t_arc(a.first, a.second)));
				}

				// a.first -> a.second -> u
				set<t_vertex> succ = set_substract(get_successors(a.second), i->second);
				succ.erase(a.first);
				for(set<t_vertex>::const_iterator u = succ.begin(); u != succ.end(); u++){
					decrease_diamond(t_diamond(a.first, *u));
					p3_set.erase(t_3path(t_arc(a.first, a.second),t_arc(a.second, *u)));
				}

				// second, add new P3 and diamonds
				set<t_vertex> belt = set_intersect(i->second,j->second);
				for(set<t_vertex>::const_iterator u = belt.begin(); u != belt.end(); u++)
					p3_set.insert(t_3path(t_arc(a.first, *u),t_arc(*u, a.second)));
				if(belt.size()>1)
					diamond_set.insert(
							pair<t_diamond,uint>(t_diamond(a.first,a.second),belt.size())
						);
				return true;
			} else return false;
		else return false;
	else return false;
}
// mark an arc forbidden/permanent,
// note that, if adjacent arcs are also forbidden/permanent, this mark may cascade
bool t_dir_graph::mark_arc(const t_arc& a, const uint mark){
	if(successors.find(a.first) == successors.end()) return false; // dont mark edges that are not in VxV
	if(successors.find(a.second) == successors.end()) return false; // dont mark edges that are not in VxV
	dbgcout << "marking " << a << (mark==MRK_NONE?"none":(mark==MRK_FORBIDDEN?"forbidden":"permanent")) << "\n";
	bool result = true;
	switch(mark){
		case MRK_NONE:
			return (forbidden_arcs.erase(a) || permanent_arcs.erase(a));
			break;
		case MRK_FORBIDDEN:
			if((!contains_arc(a)) && (forbidden_arcs.find(a) == forbidden_arcs.end())
					&& forbidden_arcs.insert(a).second) {
				set<t_vertex> belt = set_intersect(get_successors(a.first), get_predecessors(a.second));
				for(set<t_vertex>::iterator u = belt.begin(); u != belt.end(); u++){
					if(permanent_arcs.find(t_arc(a.first, *u)) != permanent_arcs.end())
						result &= mark_arc(t_arc(*u, a.second), MRK_FORBIDDEN);
					if(permanent_arcs.find(t_arc(*u, a.second)) != permanent_arcs.end())
						result &= mark_arc(t_arc(a.first, *u), MRK_FORBIDDEN);
				}
				return result;
			} else return false;
			break;
		case MRK_PERMANENT:
			if(contains_arc(a) && (permanent_arcs.find(a) == permanent_arcs.end())
					&& permanent_arcs.insert(a).second) {
				set<t_vertex> pred = get_predecessors(a.first);
				for(set<t_vertex>::iterator u = pred.begin(); u != pred.end(); u++){
					if(permanent_arcs.find(t_arc(*u, a.first)) != permanent_arcs.end())
						result &= mark_arc(t_arc(*u, a.second), MRK_PERMANENT);
				}
				pred = set_substract(get_vertices(), get_predecessors(a.second));
				for(set<t_vertex>::iterator u = pred.begin(); u != pred.end(); u++){
					if(forbidden_arcs.find(t_arc(*u, a.second)) != forbidden_arcs.end())
						result &= mark_arc(t_arc(*u, a.second), MRK_FORBIDDEN);
				}
				set<t_vertex> succ = get_successors(a.second);
				for(set<t_vertex>::iterator u = succ.begin(); u != succ.end(); u++){
					if(permanent_arcs.find(t_arc(a.second, *u)) != permanent_arcs.end())
						result &= mark_arc(t_arc(a.first, *u), MRK_PERMANENT);
				}
				succ = set_substract(get_vertices(), get_successors(a.first));
				for(set<t_vertex>::iterator u = succ.begin(); u != succ.end(); u++){
					if(forbidden_arcs.find(t_arc(a.first, *u)) != forbidden_arcs.end())
						result &= mark_arc(t_arc(a.second, *u), MRK_FORBIDDEN);
				}
				return result;
			} else return false;
			break;
		default: return false;
	}	
}
// delete all vertices that are not in the given vertex set from the graph
void t_dir_graph::induced_subgraph(const set<t_vertex>& vertices){
	set<t_vertex> v = set_substract(get_vertices(), vertices);
	for(set<t_vertex>::const_iterator i = v.begin(); i != v.end(); i++)
		delete_vertex(*i);
}

// split a weakly connected component from the digraph
// reachability by BFS
t_dir_graph* t_dir_graph::split_component(){
	set<t_vertex> comp = get_reachable(successors.begin()->first, true);
	
	t_dir_graph* result = get_induced_subgraph(comp);
	for(set<t_vertex>::const_iterator i = comp.begin(); i != comp.end(); i++)
		delete_vertex(*i);

	return result;
}

// ============= IO ====================

// code snatched from http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
void tokenize(const string& str, vector<string>& tokens, const string& delimiters = " "){
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);
	while (string::npos != pos || string::npos != lastPos){
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

bool t_dir_graph::read_format2(istream& input, const string first_line){
	vector<t_vertex> verts;
	int adj;
	t_vertex v;
	uint datas = atoi(first_line.c_str());

	while(input.peek() == '\n') input.ignore(1); // read the '\n'
	// read vertices
	for(uint i = 0; i < datas; i++){
		input >> v;
		verts.push_back(v);
		while(input.peek() == '\n') input.ignore(1); // read the '\n'
	}
	for(uint i=0; i < verts.size(); i++) insert_vertex(verts[i]);
	// read adjacency matrix
	for(uint i=0; i < verts.size(); i++){
		for(uint j=0; j < verts.size(); j++){
			input >> adj;
			if(adj > 0)
				insert_arc(t_arc(verts[i], verts[j]));
			while(input.peek() == ' ') input.ignore(1); // read the ' '
			while(input.peek() == '\n') input.ignore(1); // read the '\n'
		}
	}
	return true;
}

// read a graph from the given file, return success
bool t_dir_graph::read_format1(istream& input, const string first_line){
	string arcs;
	vector<t_vertex> vert_vec;
	vector<t_vertex> arc_vec;

	input >> arcs;
	tokenize(first_line, vert_vec, ",");
	tokenize(arcs, arc_vec, "),(");

	for(uint i=0; i < vert_vec.size(); i++)
		insert_vertex(vert_vec[i]);
	size_t arrow_pos = 0;
	for(uint i=0; i < arc_vec.size(); i++){
		if((arrow_pos = arc_vec[i].find("->")) != string::npos){
			dbgcout << "arc: " << arc_vec[i].substr(0,arrow_pos) << "," << arc_vec[i].substr(arrow_pos+2) << "\n";
			insert_arc(t_arc(arc_vec[i].substr(0,arrow_pos),arc_vec[i].substr(arrow_pos+2)));
		} else {
			i++;
			if(i < arc_vec.size()){
				dbgcout << "arc: " << arc_vec[i-1] << "," << arc_vec[i] << "\n";
				insert_arc(t_arc(arc_vec[i-1],arc_vec[i]));
			}
		}
	}

	return true;
}

bool t_dir_graph::read_from_stream(istream& input){
	string first_line;

	successors.clear();
	predecessors.clear();
	permanent_arcs.clear();
	forbidden_arcs.clear();
	p3_set.clear();
	diamond_set.clear();


	input >> first_line;
	if(first_line.find(",") != string::npos)
		return read_format1(input, first_line);
	else
		return read_format2(input, first_line);
}

// read a graph from the given file, return success
bool t_dir_graph::read_from_file(const char* filename){
	
	if(!strcmp(filename,"stdin")) read_from_stream(cin); else {
		ifstream* file = new ifstream(filename, ios::in|ios::binary);
		
		if(!file) {
			cout << "error opening \"" << filename << "\"\n";
			return false;
		}
		dbgcout << "reading data from " << filename << "\n";
		read_from_stream(*file);
		delete file;
	}
	return true;
}

// print the digraph to a stream
void t_dir_graph::print(ostream& output) const{
	set<t_vertex> v;
	set<t_arc> a;
	for(t_adjtable::const_iterator i = successors.begin(); i != successors.end(); i++){
		v.insert(i->first);
		for(t_adjlist::const_iterator j = i->second.begin(); j != i->second.end(); j++){
			a.insert(t_arc(i->first,*j));
		}
	}
	// print vertex set
	output << "(";
	for(set<t_vertex>::const_iterator i = v.begin(); i != v.end(); i++)
		output << *i << ",";
	output << "\b)\n";
	// print arcs
	for(set<t_arc>::const_iterator i = a.begin(); i != a.end(); i++)
		output << *i << ((permanent_arcs.find(*i) == permanent_arcs.end())?",":"!,");
	output << "\b";
}
// streaming IO
ostream& operator<<(ostream& os, const t_dir_graph& d){
	d.print(os);
	return os;
}

istream& operator>>(istream& is, t_dir_graph& d){
	d.read_from_stream(is);
	return is;
}

// ========== end Class ============

#define ARC_DELIMETER "->"
ostream& operator<<(ostream& os, const t_arc& a){
	os << "(" << a.first << ARC_DELIMETER << a.second << ")";
	return os;
}

t_dir_graph* generate_random_digraph(const uint vertices, const double arc_probability){
	t_dir_graph* d = new t_dir_graph();
	srand(time(NULL));
	char vname[100];
	char vname2[100];
	for(uint i = 0; i < vertices; i++) {
		sprintf( vname, "%d", i );
		d->insert_vertex(vname);
	}
	for(uint i = 0; i < vertices; i++)
		for(uint j = 0; j < vertices; j++) if(i!=j){
			sprintf( vname, "%d", i );
			sprintf( vname2, "%d", j );
			if(rand() < arc_probability * RAND_MAX) d->insert_arc(t_arc(vname,vname2));
		}
	return d;
}


#endif

/***************************************************
 * editing.cpp
  * by M.Weller
 **************************************************/
#ifndef editing_cpp
#define editing_cpp

#include <ctime>
#include "editing.h"
#include <iostream>


ostream& operator<<(ostream& os, const t_3path& p){
	os << "(" << p.first.first << "->" << p.first.second << "->" << p.second.second << ")";
	return os;
}

ostream& operator<<(ostream& os, const set<t_3path>& sp){
	if(sp.size()){
		os << "(";
		for(set<t_3path>::const_iterator i = sp.begin(); i != sp.end(); i++)
			os << *i << ",";
		os << "\b)";
	} else os << "()";
	return os;
}

ostream& operator<<(ostream& os, const t_edit& e){
	os << (e.first?"+":"-") << e.second;
	return os;
}

ostream& operator<<(ostream& os, const t_edit_ex& e_ex){
	t_edit e = e_ex.first;
	os << e;
	if(e_ex.second.length()) os << "[" << e_ex.second << "]";
	return os;
}
ostream& operator<<(ostream& os, const t_edit_ex_history& h){
	for(t_edit_ex_history::const_iterator i = h.begin(); i != h.end(); i++)
		os << *i << ",";
	os << "\b";
	return os;
}

void append_history(t_edit_ex_history* h1, const t_edit_ex_history* h2){
	if(h1)
		if(h2)
			for(t_edit_ex_history::const_iterator i = h2->begin(); i != h2->end(); i++)
				h1->push_back(*i);
}

inline bool p3_disjoint(const t_3path& p, const t_3path& q, const bool allow_insert){
	if(allow_insert){
		t_arc zp,zq;
		zp = t_arc(p.first.first, p.second.second);
		zq = t_arc(q.first.first, q.second.second);
		return (p.first != q.first) && (p.first != q.second) &&
				(p.second != q.first) && (p.second != q.second) &&
				(zp != zq);
	} else return (p.first != q.first) && (p.first != q.second) &&
				(p.second != q.first) && (p.second != q.second);
}

#include <algorithm>
// sort a (p3,(#cliques,|cliques|)) pair #cliques-biased
bool my_cmp(const pair<t_3path, pair<uint, uint> >& p1,
			const pair<t_3path, pair<uint, uint> >& p2){
	if(p1.second.first < p2.second.first) return true;
	if(p1.second.first > p2.second.first) return false;
	// from here, we know that p1.second.first == p2.second.first
	if(p1.second.second < p2.second.second) return true;
		else return false;
}
// return the P3 with the highest clique-number, that is, the one conflicting with
// the most other P3
p_3path get_branching_P3(const t_dir_graph& d, const set<t_3path>& paths, const bool allow_insert){
	map<t_arc, set<t_3path> > conf_cliques;
	map<t_3path, pair<uint,uint> > clique_count;
	set<t_arc> forbidden_arcs;

	// 1. fill the data structures
	set<t_vertex> vertices = d.get_vertices();
	for(set<t_vertex>::const_iterator u = vertices.begin(); u != vertices.end(); u++)
		for(set<t_vertex>::const_iterator v = vertices.begin(); v != vertices.end(); v++) if(*u != *v)
			conf_cliques.insert(pair<t_arc,set<t_3path> >(t_arc(*u,*v),set<t_3path>()));
	// each arc induces a clique in the conflict graph
	for(set<t_3path>::const_iterator p = paths.begin(); p != paths.end(); p++){
		clique_count.insert(pair<t_3path, pair<uint,uint> >(*p,pair<uint,uint>(0,0)));
		conf_cliques[p->first].insert(*p);
		conf_cliques[p->second].insert(*p);
		conf_cliques[t_arc(p->first.first,p->second.second)].insert(*p);
	}
	// count cliques that are of size > 1 and their sizes
	for(map<t_arc, set<t_3path> >::const_iterator i = conf_cliques.begin(); i != conf_cliques.end(); i++)
		if(i->second.size() > 1) 
			for(set<t_3path>::const_iterator p = i->second.begin(); p != i->second.end(); p++){
				clique_count[*p].first++;
				clique_count[*p].second += i->second.size();
			}
	p_3path result = NULL;
	pair<uint, uint> max_val = pair<uint,uint>(0,0);
	for(map<t_3path, pair<uint,uint> >::const_iterator i = clique_count.begin(); i != clique_count.end(); i++)
		if(result){
			if(my_cmp(pair<t_3path,pair<uint,uint> >(*result, max_val),*i)) max_val = i->second;
		} else result = new t_3path(i->first);
	return result;
}
// O(n^3 log n)
set<t_3path> get_maximal_disjoint_P3(const t_dir_graph& d, const set<t_3path>& paths, const bool allow_insert){
	map<t_arc, set<t_3path> > conf_cliques;
	map<t_3path, pair<uint,uint> > clique_count;
	set<t_arc> forbidden_arcs;

	// 1. fill the data structures
	set<t_vertex> vertices = d.get_vertices();
	for(set<t_vertex>::const_iterator u = vertices.begin(); u != vertices.end(); u++)
		for(set<t_vertex>::const_iterator v = vertices.begin(); v != vertices.end(); v++) if(*u != *v)
			conf_cliques.insert(pair<t_arc,set<t_3path> >(t_arc(*u,*v),set<t_3path>()));
	// each arc induces a clique in the conflict graph
	for(set<t_3path>::const_iterator p = paths.begin(); p != paths.end(); p++){
		clique_count.insert(pair<t_3path, pair<uint,uint> >(*p,pair<uint,uint>(0,0)));
		conf_cliques[p->first].insert(*p);
		conf_cliques[p->second].insert(*p);
		conf_cliques[t_arc(p->first.first,p->second.second)].insert(*p);
	}
	// count cliques that are of size > 1 and their sizes
	for(map<t_arc, set<t_3path> >::const_iterator i = conf_cliques.begin(); i != conf_cliques.end(); i++)
		if(i->second.size() > 1) 
			for(set<t_3path>::const_iterator p = i->second.begin(); p != i->second.end(); p++){
				clique_count[*p].first++;
				clique_count[*p].second += i->second.size();
			}
	
	// 2. sort by clique count & size
	vector<pair<t_3path, pair<uint, uint> > > count_vec(clique_count.begin(), clique_count.end());
	sort(count_vec.begin(), count_vec.end(), my_cmp);
	set<t_3path> result;
	// 3. one by one, take the minimum elements:
	// 	if they are forbidden, discard them,
	// 	else insert into the result and mark their 3 arcs forbidden for the future
	for(vector<pair<t_3path, pair<uint, uint> > >::iterator i = count_vec.begin(); i != count_vec.end(); i++)
		if((forbidden_arcs.find(i->first.first) == forbidden_arcs.end()) &&
			(forbidden_arcs.find(i->first.second) == forbidden_arcs.end()) &&
			(forbidden_arcs.find(t_arc(i->first.first.first,i->first.second.second)) == forbidden_arcs.end())){
			result.insert(i->first);
			forbidden_arcs.insert(i->first.first);
			forbidden_arcs.insert(i->first.second);
			forbidden_arcs.insert(t_arc(i->first.first.first,i->first.second.second));
		}
	return result;
}
// apply preprocessing rule0: process high-degree rules
uint prep_high_degree(
		t_dir_graph* d,
		const uint prev_k,
		const bool allow_insert,
		t_edit_ex_history* h){
	if(prev_k >= MY_INFINITY) return prev_k;
	set<t_vertex> vertices = d->get_vertices();
	set<t_vertex> pred, succ;
	uint k = prev_k;
	uint k_mod = 0;
	t_3path q;
	set<t_3path> p3_temp, to_delete;

	set<t_3path> disjoint_p3 = get_maximal_disjoint_P3(*d, d->get_p3(), allow_insert);
	if((int)k - (int)disjoint_p3.size() < 0) return MY_INFINITY;

	for(set<t_vertex>::const_iterator i = vertices.begin(); i != vertices.end(); i++)
		for(set<t_vertex>::const_iterator j = vertices.begin(); j != vertices.end(); j++) if(i != j){
			// modify k by the number of disjoint p3
			k_mod = disjoint_p3.size() - 1;
			if((!d->contains_arc(t_arc(*i,*j)))){
				if(allow_insert){
					// a)  more then k common neighbors: insert arc
					succ = d->get_successors(*i);
					pred = d->get_predecessors(*j);
					p3_temp = disjoint_p3;
//					cout << t_arc(*i,*j) << ": "<< pred << ", "<<succ<<" disjoint_p3: "<<disjoint_p3<<"\n";
					to_delete.clear();
					for(set<t_vertex>::const_iterator u = pred.begin(); u != pred.end(); u++)
						for(set<t_3path>::const_iterator p = p3_temp.begin(); p != p3_temp.end(); p++)
							if(!p3_disjoint(*p,t_3path(t_arc(*i,*u),t_arc(*u,*j)),allow_insert))
								to_delete.insert(*p);
					for(set<t_3path>::const_iterator p = to_delete.begin(); p != to_delete.end(); p++)
						p3_temp.erase(*p);
			
					to_delete.clear();
//					cout << t_arc(*i,*j) << ": "<< pred << ", "<<succ<<"\n";
					for(set<t_vertex>::const_iterator u = succ.begin(); u != succ.end(); u++)
						for(set<t_3path>::const_iterator p = disjoint_p3.begin(); p != disjoint_p3.end(); p++)
							if(!p3_disjoint(*p,t_3path(t_arc(*i,*u),t_arc(*u,*j)),allow_insert))
								to_delete.insert(*p);
					for(set<t_3path>::const_iterator p = to_delete.begin(); p != to_delete.end(); p++)
						p3_temp.erase(*p);
					k_mod = p3_temp.size();
	
					if(set_intersect(pred, succ).size()+k_mod > k){
						if(k > k_mod){
							if(d->insert_arc(t_arc(*i,*j))) {
								k--;
								if(h) h->push_back(t_edit_ex(t_edit(ED_INSERT, t_arc(*i,*j)),"H"));
								prepdbgcout << "high degree: " << t_arc(*i,*j) << " is inserted, new k: "<<k<<"\n";
								if(k==0) return 0;
							}
							disjoint_p3 = get_maximal_disjoint_P3(*d, d->get_p3(), allow_insert);
						} else {
							prepdbgcout << "high degree: could not insert " << t_arc(*i,*j) << "\n";
							return MY_INFINITY;
						}
					}
				}
			} else {
				// b)  more then k uncommon neighbors: delete arc
				// only consider vertices that are part of a P3 containing (*i,*j)
				pred = set_substract(d->get_predecessors(*i), d->get_predecessors(*j));
				succ = set_substract(d->get_successors(*j), d->get_successors(*i));
				
				pred.erase(*j);
				succ.erase(*i);

				p3_temp = disjoint_p3;
//				cout << t_arc(*i,*j) << ": "<< pred << ", "<<succ<<" disjoint_p3: "<<disjoint_p3<<"\n";
				to_delete.clear();
				for(set<t_vertex>::const_iterator u = pred.begin(); u != pred.end(); u++)
					for(set<t_3path>::const_iterator p = p3_temp.begin(); p != p3_temp.end(); p++)
						if(!p3_disjoint(*p,t_3path(t_arc(*u,*i),t_arc(*i,*j)),allow_insert))
							to_delete.insert(*p);
				for(set<t_3path>::const_iterator p = to_delete.begin(); p != to_delete.end(); p++)
					p3_temp.erase(*p);
	
				to_delete.clear();
				for(set<t_vertex>::const_iterator u = succ.begin(); u != succ.end(); u++)
					for(set<t_3path>::const_iterator p = disjoint_p3.begin(); p != disjoint_p3.end(); p++)
						if(!p3_disjoint(*p,t_3path(t_arc(*i,*j),t_arc(*j,*u)),allow_insert))
							to_delete.insert(*p);
				for(set<t_3path>::const_iterator p = to_delete.begin(); p != to_delete.end(); p++)
					p3_temp.erase(*p);

//				cout << t_arc(*i,*j) << ":  disjoint p3: " << p3_temp << "\n";
				k_mod = p3_temp.size();

				if(set_unite(pred, succ).size() + k_mod > k){
					if(k > k_mod){
						if(d->delete_arc(t_arc(*i,*j))) {
							k--;
							if(h) h->push_back(t_edit_ex(t_edit(ED_DELETE, t_arc(*i,*j)),"H"));
							prepdbgcout << "high degree: " << t_arc(*i,*j) << " is deleted, new k: "<<k<<"\n";
							if(k==0) return 0;
						}
						disjoint_p3 = get_maximal_disjoint_P3(*d, d->get_p3(), allow_insert);
					} else {
						prepdbgcout << "high degree: could not delete " << t_arc(*i,*j) << "\n";
						return MY_INFINITY;
					}
				}
			}
	}
	return k;
}

// apply preprocessing rule1: remove all good vertices
void prep_remove_good(t_dir_graph* d){
	set<t_vertex> vertices(d->get_vertices());
	set<t_vertex> u_adjacent;
	set<t_vertex> v_adjacent;
	set<t_vertex> dirty;
	
	for(set<t_vertex>::const_iterator u = vertices.begin();u != vertices.end();u++){
		u_adjacent = d->get_successors(*u);
		for(set<t_vertex>::const_iterator v = u_adjacent.begin();v != u_adjacent.end();v++){
			v_adjacent = d->get_successors(*v);
			for(set<t_vertex>::const_iterator w = v_adjacent.begin();w != v_adjacent.end();w++){
				if(is_3path(*d,*u,*v,*w)) {
					dirty.insert(*u);
					dirty.insert(*v);
					dirty.insert(*w);
				}
			}
		}
	}
	d->induced_subgraph(dirty);
	prepdbgcout << "after removing clean: " << *d << "\n";
}

// functions for splitting and rejoining vertices that have a single double-arc
// and indeg = 1 or outdeg = 1, that is, no other incoming or no other outgoing arcs
void split_double_end(t_dir_graph* d, t_vertex end){
	set<t_vertex> pred = d->get_predecessors(end);
	set<t_vertex> succ = d->get_successors(end);
	set<t_vertex> inter = set_intersect(pred, succ);

	if((inter.size()==1) && ((succ.size()==1) || (pred.size()==1))){
		t_vertex new_end = end + SPLIT_SEPERATOR;
		d->insert_vertex(new_end);
		d->insert_arc(t_arc(new_end, end));

		d->delete_arc(t_arc(new_end, *(inter.begin())));
		if(pred.size()==1){
			succ.erase(*(inter.begin()));
			for(set<t_vertex>::iterator i = succ.begin(); i != succ.end(); i++)
				d->insert_arc(t_arc(new_end, *i));
		} else {
			pred.erase(*(inter.begin()));
			for(set<t_vertex>::iterator i = pred.begin(); i != pred.end(); i++)
				d->insert_arc(t_arc(*i, new_end));
		}
	}
}

void join_split_vertices(t_dir_graph* d){
	set<t_vertex> vertices = d->get_vertices();
	set<t_vertex> succ, pred;
	set<t_vertex> to_remove;
	t_vertex end1, end2;
	size_t pos;
	for(set<t_vertex>::iterator i = vertices.begin(); i != vertices.end(); i++)
		if((pos = i->find(SPLIT_SEPERATOR,0)) != string::npos){
			end1 = *i;
			end1.erase(pos + 1);
			end2 = end1;
			end2.erase(pos + 1 - sizeof(SPLIT_SEPERATOR));

			to_remove.insert(end1);
			pred = d->get_predecessors(end1);
			succ = d->get_successors(end1);
			for(set<t_vertex>::iterator j = pred.begin(); j != pred.end(); j++)
				if(*j != end2) d->insert_arc(t_arc(*j, end2));
			for(set<t_vertex>::iterator j = succ.begin(); j != succ.end(); j++)
				if(*j != end2) d->insert_arc(t_arc(end2, *j));
		}
	for(set<t_vertex>::iterator i = to_remove.begin(); i != to_remove.end(); i++)
		d->delete_vertex(*i);
}

void remove_double_arc(t_dir_graph* d, t_vertex v1, t_vertex v2){
	t_vertex u,v;
	set<t_vertex> pre_u = d->get_predecessors(v1);
	set<t_vertex> suc_u = d->get_successors(v1);
	set<t_vertex> pre_v = d->get_predecessors(v2);
	set<t_vertex> suc_v = d->get_successors(v2);
	set<t_vertex> temp;
	if((pre_u.find(v2) != pre_u.end()) && (suc_u.find(v2) != suc_u.end())){
		if((pre_u.size() == 1) && (suc_v.size() == 1)){
			u = v1;
			v = v2;
		} else if((pre_v.size() == 1) && (suc_u.size() == 1)){
			v = v1;
			u = v2;
			// swap pre_u/v and suc_u/v respectively
			temp = pre_u;pre_u = pre_v;pre_v = temp;
			temp = suc_u;suc_u = suc_v;suc_v = temp;
		} else return;
		// from here, we know that indeg(u)==1 and outdeg(v)==1
		// if removing (u,v) does not create any P3, remove it from consideration
		if(!(set_intersect(pre_v, suc_u).size()))
			d->delete_arc(t_arc(u,v));
	}
}

// apply preprocessing rule: remove doublearc-parts where possible
void prep_remove_double(t_dir_graph* d){
	set<t_arc> arcs = d->get_arcs();
	set<t_arc>::iterator i = arcs.begin();
	// two possibilities to remove double-arcs:
	while(i != arcs.end()) if(arcs.find(reverse_arc(*i)) != i){
//		cout << "checking arc " << *i << "...\n";
		// first: remove a part of the double-arc where we can
		remove_double_arc(d, i->first, i->second);
		// second: at the end of a graph, split the terminating vertex in a source and a sink
		split_double_end(d, i->first);
		split_double_end(d, i->second);
		i++;
	}
}



#define DIA_SINKS true
#define DIA_SOURCES false

bool sources_diamond_check(const t_dir_graph* d, const set<t_vertex> diamond_middle, const set<t_vertex> sources){
	uint single_succ = 0;
	uint multi_pred = 0;
	for(set<t_vertex>::iterator m = diamond_middle.begin(); m != diamond_middle.end(); m++){
		if(set_intersect(d->get_predecessors(*m), sources).size() != 1) multi_pred++;
		if(d->get_successors(*m).size() == 1) single_succ++;
	}
	return ((single_succ) && (!multi_pred)); 
	// at least one vertex in diamond_middle has a single successor
	// and none have multiple predecessors
}
bool sinks_diamond_check(const t_dir_graph* d, const set<t_vertex> diamond_middle, const set<t_vertex> sinks){
	uint single_pred = 0;
	uint multi_succ = 0;
	for(set<t_vertex>::iterator m = diamond_middle.begin(); m != diamond_middle.end(); m++){
		if(d->get_predecessors(*m).size() == 1) single_pred++;
		if(set_intersect(d->get_successors(*m), sinks).size() != 1) multi_succ++;
	}
	return ((single_pred) && (!multi_succ));
	// at least one vertex in diamond_middle has a single predecessor
	// and none have multiple successors
}
// detects an illegal structure: an n-diamond with ( <n )-sided tentacles:
// e.g. 3-diamond with 2-sided tentacles:
/*     o
//    /|\
//   o o o
//  /|\|/|\
// o o o o o
*/
set<t_vertex> remove_illegal(t_dir_graph* d,set<t_vertex>* pred, set<t_vertex>* N, set<t_vertex>* succ, bool sink){
	set<t_vertex> removal;
	set<t_vertex> common_pred, common_succ, succ1, succ2, pred1, pred2;

	for(set<t_vertex>::iterator i = N->begin(); i != N->end(); i++)
		for(set<t_vertex>::iterator j = i; j != N->end(); j++) if(i != j){
			succ1 = d->get_successors(*i);
			pred1 = d->get_predecessors(*i);
			succ2 = d->get_successors(*j);
			pred2 = d->get_predecessors(*j);

			common_pred.clear(); common_succ.clear();
			common_pred = set_intersect(pred1, pred2);
			common_succ = set_intersect(succ1, succ2);

			if(pred) common_pred = set_intersect(common_pred, *pred); // limit to the given sources
			if(succ) common_succ = set_intersect(common_succ, *succ); // limit to the given sinks

			if(common_pred.size() && common_succ.size()){
//				t_vertex diamond_pred = *(common_pred.begin());
//				t_vertex diamond_succ = *(common_succ.begin());
//				set<t_vertex> diamond_middle;
//				diamond_middle = set_intersect(d->get_successors(diamond_pred), d->get_predecessors(diamond_succ));

				eddbgcout << "diamond: " << common_pred << ", (" << *i << ", " << *j << "), " << common_succ << "\n";
				if(sink == DIA_SINKS){
				//	if(sinks_diamond_check(d, diamond_middle, (succ?*succ:set<t_vertex>::set<t_vertex>())))
						// found an illegal structure, mark for removal from N
						removal.insert(*i);
						removal.insert(*j);
//						removal = set_unite(removal, diamond_middle);
				} else {
				//	if(sources_diamond_check(d, diamond_middle, (pred?*pred:set<t_vertex>::set<t_vertex>()))) 
						// found an illegal structure, mark for removal from N
						removal.insert(*i);
						removal.insert(*j);
//						removal = set_unite(removal, diamond_middle);
				}
				eddbgcout << "for removal now: " << removal << "\n";
			}
		}
	return removal;
}

#include <sstream>

// apply preprocessing rule2: cut the endings and return the new k
uint prep_cut_endings(
		t_dir_graph* d,
		const uint prev_k,
		const bool allow_insert,
		t_edit_ex_history* h){
	if(prev_k >= MY_INFINITY) return prev_k;
	// remember to remove isolated diamonds from the structure !
	set<t_vertex> N, pred, succ;
	set<t_vertex> diamond_infested;
//	set<t_vertex> sinks(d->get_sinks());
//	set<t_vertex> sources(d->get_sources());
//	set<t_vertex> sources = set_substract(d->get_p3_first(),set_unite(d->get_p3_second(), d->get_p3_third()));
//	set<t_vertex> sinks = set_substract(d->get_p3_third(),set_unite(d->get_p3_second(), d->get_p3_first()));
	set<t_vertex> s, first, second, third;

	t_dir_graph* temp;
	uint k = prev_k;
	bool changed = true;
	stringstream explain;
	set<t_vertex> d_free,T,R,R_t;
	uint maxx_R_t;

	while(changed){
		temp = new t_dir_graph(*d);
		if(allow_insert)
			s = d->get_sinks();
		else 
			s = set_substract(d->get_vertices(), d->get_p3_second());
//		first = d->get_p3_first();
//		third = d->get_p3_third();
		changed = false;
		// process sinks
		if(s.size()){
			N.clear();
			for(set<t_vertex>::const_iterator sink = s.begin(); sink != s.end(); sink++)
				N = set_unite(N, temp->get_predecessors(*sink));
			prepdbgcout << "processing sinks " << s << " and their neighbors "<<N<<"\n";
//			remove_diamonds(temp, NULL, &N, &sinks);
			if(allow_insert){
				diamond_infested = remove_illegal(temp, NULL, &N, &s, DIA_SINKS); // instead of remove_diamonds
				d_free = set_substract(N, diamond_infested);
			} else d_free = N;

			for(set<t_vertex>::const_iterator i = d_free.begin(); i != d_free.end(); i++){
				T.clear(); R.clear();
				succ = set_intersect(s, temp->get_successors(*i));
				pred = set_substract(temp->get_predecessors(*i), succ);
				
				prepdbgcout << "pred= " << pred << " succ="<<succ<<"\n";
				// for each a in R, there is a u in succ(v) with (a,u) not in A
				for(set<t_vertex>::const_iterator j = succ.begin(); j != succ.end(); j++)
					if(!(set_substract(pred, temp->get_predecessors(*j)).empty())) R.insert(*j);
				// for each t in T, there is an a in R with (a,t) not in A
				for(set<t_vertex>::const_iterator j = pred.begin(); j != pred.end(); j++)
					if(!(set_substract(R, temp->get_successors(*j)).empty())) T.insert(*j);
				if(T.size() && (T.size() <= R.size())){
					maxx_R_t = 0;
					set<t_vertex> R_t;
					// only apply the rule, if forall T`<=T |R|>=|T`|+|bigcap_{t\in T`} R_t|
					// since this is too hard to check. Hence, only apply, iff
					// |R| >= |T| + max_{t\in T}|R_t|
					for(set<t_vertex>::const_iterator j = T.begin(); j != T.end(); j++){
						R_t = set_intersect(temp->get_successors(*j),R);
						if(R_t.size() > maxx_R_t) maxx_R_t = R_t.size();
					}
	
					if((T.size() + maxx_R_t <= R.size())){
						prepdbgcout << "R= " << R << " T="<<T<<"\n";
						for(set<t_vertex>::const_iterator j = T.begin(); j != T.end(); j++)
							if(k) {
								if(d->delete_arc(t_arc(*j,*i), true)) {
									changed = true;
									k--;
									explain.str("");
									explain << "Si(" << *i << ")";
									if(h) h->push_back(t_edit_ex(t_edit(ED_DELETE, t_arc(*j,*i)),explain.str()));
									prepdbgcout << "sinks: " << t_arc(*j,*i) << " is deleted\n";
									if(k==0) return 0;
								}
							} else {
								prepdbgcout << "sinks: could not delete " << t_arc(*j,*i) << "\n";
								delete temp;
								return MY_INFINITY;
							}
					}
				}
			}
		}
		if(!changed){
			if(allow_insert)
				s = d->get_sources();
			if(s.size()){
				delete temp;temp = new t_dir_graph(*d); N.clear();
				// process sources
				prepdbgcout << "processing sources " << s << "\n";
				N.clear();
				for(set<t_vertex>::const_iterator source = s.begin(); source != s.end(); source++)
					N = set_unite(N, temp->get_successors(*source));
	//			remove_diamonds(temp, &sources, &N, NULL);
				if(allow_insert){
					diamond_infested = remove_illegal(temp, &s, &N, NULL, DIA_SOURCES); // instead of remove_diamonds
					d_free = set_substract(N, diamond_infested);
				} else d_free = N;
				
				for(set<t_vertex>::const_iterator i = d_free.begin(); i != d_free.end(); i++){
					T.clear(); R.clear();
					pred = set_intersect(s, temp->get_predecessors(*i));
					succ = set_substract(temp->get_successors(*i), pred);
	
					// for each a in R, there is a u in succ(v) with (a,u) not in A
					for(set<t_vertex>::const_iterator j = pred.begin(); j != pred.end(); j++)
						if(!(set_substract(succ, temp->get_successors(*j)).empty())) R.insert(*j);
					// for each t in T, there is an a in R with (a,t) not in A
					for(set<t_vertex>::const_iterator j = succ.begin(); j != succ.end(); j++)
						if(!(set_substract(R, temp->get_predecessors(*j)).empty())) T.insert(*j);
					if(T.size() && (T.size() <= R.size())){
						maxx_R_t = 0;
						// only apply the rule, if forall T`<=T |R|>=|T`|+|bigcap_{t\in T`} R_t|
						// since this is too hard to check, only apply, if
						// |R| >= |T| + max_{t\in T}|R_t|
						for(set<t_vertex>::const_iterator j = T.begin(); j != T.end(); j++){
							R_t = set_intersect(temp->get_predecessors(*j),R);
							if(R_t.size() > maxx_R_t) maxx_R_t = R_t.size();
						}
						if((T.size() + maxx_R_t <= R.size())){
							for(set<t_vertex>::const_iterator j = T.begin(); j != T.end(); j++)
								if(k){
									if(d->delete_arc(t_arc(*i,*j), true)) {
										changed = true;
										k--;
										explain.str("");
										explain << "So(" << *i << ")";
										if(h) h->push_back(t_edit_ex(t_edit(ED_DELETE, t_arc(*i,*j)),explain.str()));
										prepdbgcout << "sources: " << t_arc(*i,*j) << " is deleted\n";
										if(k==0) return 0;
									}
								} else {
									prepdbgcout << "sources: could not delete " << t_arc(*i,*j) << "\n";
									delete temp;
									return MY_INFINITY;
								}
						}
					}
				}
			}
		}
		delete temp;
	}
	return k;
}

// solve transitive vertx delete problem
uint trans_vertex_del(t_dir_graph* d){
	set<t_vertex> v = d->get_vertices();
	if(!get_3paths(*d).empty()){
		t_dir_graph* temp;
		t_dir_graph* best = NULL;
		uint k = MY_INFINITY;
		uint newk = k;
		for(set<t_vertex>::iterator i = v.begin(); i != v.end(); i++){
			// make a copy of the graph
			temp = new t_dir_graph(*d);
			// remove the current vertex
			temp->delete_vertex(*i);
			// branch
			newk = trans_vertex_del(temp);
			if(newk < k){
				k = newk;
				if(best) delete best;
				best = new t_dir_graph(*temp);
			}
			delete temp;
		}
		*d = *best;
		delete best;
		return 1 + k;
	} else return 0;
}
// filter arcs beginning in 'start' and ending in 'end' from arc set A
// if 'start' or 'end' is NULL, allow all vertices as start or endpoints respectively
set<t_arc> arc_filter(const set<t_arc> A, const set<t_vertex> start, const set<t_vertex> end){
	set<t_arc> result;
	for(set<t_arc>::iterator i = A.begin(); i != A.end(); i++)
		if(start.find(i->first)!=start.end())
			if(end.find(i->second)!=end.end())
				result.insert(*i);
	return result;
}

//  apply all preprocessing rules with k = prev_k, meaning the question is:
//  can d be turned transitive with at most prev_k edits ?
uint apply_preprocessing(
		t_dir_graph* d,
		const uint prev_k,
		const bool allow_insert,
		t_edit_ex_history* h){
	uint corr_k = prev_k;
	uint tmp;
	do{
		do{
			do{
				tmp = corr_k;
				prepdbgcout << "=== applying prep-rule1: remove_good (k = " << corr_k << ")===\n";
//				join_split_vertices(d);
				prep_remove_good(d);
				prepdbgcout << "=== applying prep-rule1a: remove_double_arcs (k = " << corr_k << ")===\n";
//				prep_remove_double(d);
				prepdbgcout << "=== applying prep-rule2: high_degree (k = " << corr_k << ")===\n";
				corr_k = prep_high_degree(d, corr_k, allow_insert, h);
			} while((corr_k < MY_INFINITY) && (tmp != corr_k));
			tmp = corr_k;
			prepdbgcout << "=== applying prep-rule3: cut_endings (k = " << corr_k << ")===\n";
			corr_k = prep_cut_endings(d, corr_k, allow_insert, h);
		} while((corr_k < MY_INFINITY) && (tmp != corr_k));
		tmp = corr_k;
//		prepdbgcout << "=== applying prep-rule4: crown_reduct(k = " << corr_k << ")===\n";
//		corr_k = prep_crown_reduct(d, corr_k, h);
	} while((corr_k < MY_INFINITY) && (tmp != corr_k));
	// undo the arc split done by prep_remove_double
	join_split_vertices(d);
	return corr_k;
}


bool is_3path(const t_dir_graph& d,const t_vertex& u,const t_vertex& v,const t_vertex& w){
	if((u != v) && (u != w) && (v != w)) {
		set<t_vertex> adj = d.get_successors(u);
		if((adj.find(v) != adj.end()) && (adj.find(w) == adj.end())) {
			adj = d.get_successors(v);
			if(adj.find(w) != adj.end()) return true;
				else return false;
		} else return false;
	} else return false;
}

// return a set of induced P3 that is in d containing at max max_p3 P3
set<t_3path> get_3paths(const t_dir_graph& d, const uint max_p3){
	return get_diamond(d, max_p3, 1);
}

// returns a set of (>min_n-1)-diamonds represented as one of its 3pathes in d of size at most max_diamonds
set<t_3path> get_diamond(const t_dir_graph& d, uint max_diamonds, uint min_n){
	set<t_3path> result;

	set<t_vertex> vertices = d.get_vertices();
	set<t_vertex>::const_iterator u = vertices.begin();
	
	set<t_vertex> u_adjacent;
	set<t_vertex>::const_iterator v;
	
	set<t_vertex> v_adjacent;
	set<t_vertex>::const_iterator w;
	
	while((u != vertices.end()) && (result.size() < max_diamonds)){
		u_adjacent = d.get_successors(*u);
		v = u_adjacent.begin();
		while((v != u_adjacent.end()) && (result.size() < max_diamonds)){
			v_adjacent = d.get_successors(*v);
			w = v_adjacent.begin();
			while((w != v_adjacent.end()) && (result.size() < max_diamonds)){
				if(is_3path(d,*u,*v,*w)) {
					if(set_intersect(d.get_predecessors(*w),d.get_successors(*u)).size() >= min_n)
						result.insert(t_3path(t_arc(*u,*v),t_arc(*v,*w)));
				}
				if(result.size() < max_diamonds) w++;
			}			
			if(result.size() < max_diamonds) v++;
		}
		if(result.size() < max_diamonds) u++;
	}
	return result;
}

// try to turn d transitive with at most prev_k edits
uint get_transitive_min_edit_distance(
		t_dir_graph& d,
		const bool allow_insert,
		t_edit_ex_history* hist,
		const uint prev_k){
	eddbgcout << "=====> testing digraph (k = "<< prev_k << ") " << d << "\n";
	if(hist) eddbgcout << "k = " << prev_k << "\tedits: " << *hist << "\n";
	uint prep_edits = 0;
	uint max_k = prev_k;

	t_edit_ex_history* prep_h = NULL;
	if(hist) prep_h = new t_edit_ex_history();
/*	max_k = apply_preprocessing(&d, prev_k, prep_h);
	if(max_k >= MY_INFINITY) {
		if(hist) delete prep_h;
		return max_k;
	} // preprocessing figured that it's not possible with prev_k edits
	prep_edits = prev_k - max_k;
	if(hist) append_history(hist, prep_h);
*/	
	set<t_3path> paths = get_3paths(d);
	eddbgcout << "remaining P3: " << paths << "\n";
	if(!paths.empty()) {// if its not transitive, branch all 3 possibilities to destroy the P3
		t_3path path = *(paths.begin());
		if(max_k == 0) { // if we dont have any edits left, return infinity
			if(hist) delete prep_h;
			eddbgcout << "can not destroy P3: " << path <<"\n";
			return MY_INFINITY;
		}
		eddbgcout << "destroying P3: " << path <<"\n";
		bool success = false;
		uint branch_edits = 0;
		t_dir_graph* d_bak = new t_dir_graph(d);
		t_edit_ex_history* h_bak = NULL;
		t_edit_ex_history* best_h = NULL;
		if(hist) {
			h_bak = new t_edit_ex_history(*hist);
			best_h = new t_edit_ex_history();
		}
		// 1. insert the missing arc and branch
		if(allow_insert){
			if(d_bak->insert_arc(t_arc(path.first.first, path.second.second),true)){
				// phase1: modify the history and branch
				if(hist) h_bak->push_back(t_edit_ex(t_edit(ED_INSERT,t_arc(path.first.first, path.second.second)),""));
				branch_edits = 1 + get_transitive_min_edit_distance(*d_bak, allow_insert, h_bak, max_k - 1);
				// phase2: evaluate the result
				if(branch_edits <= max_k){
					max_k = branch_edits;
					if(hist) *best_h = *h_bak;
				}
				success |= (branch_edits <= max_k);
				// phase3: mark arc in d
				d.mark_arc(t_arc(path.first.first, path.second.second), MRK_FORBIDDEN);
				// phase4: restore copies
				*d_bak = d;
				if(hist) *h_bak = *hist;
			} else eddbgcout << t_arc(path.first.first, path.second.second) << " is forbidden\n";
		}
		// 2. delete the first arc and branch
		if(d_bak->delete_arc(path.first, true)){
			// phase1: modify the history and branch
			if(hist) h_bak->push_back(t_edit_ex(t_edit(ED_DELETE,path.first),""));
			branch_edits = 1 + get_transitive_min_edit_distance(*d_bak, allow_insert, h_bak, max_k - 1);
			// phase2: evaluate the result
			if(branch_edits <= max_k){
				max_k = branch_edits;
				if(hist) *best_h = *h_bak;
			}
			success |= (branch_edits <= max_k);
			// phase3: mark arc in d
			d.mark_arc(path.first, MRK_PERMANENT);
			// phase4: restore copies
			*d_bak = d;
			if(hist) *h_bak = *hist;
		} else eddbgcout << t_arc(path.first) << " is permanent\n";
		// 3. delete the second arc and branch
		if(d_bak->delete_arc(path.second, true)){
			// phase1: modify the history and branch
			if(hist) h_bak->push_back(t_edit_ex(t_edit(ED_DELETE,path.second),""));
			branch_edits = 1 + get_transitive_min_edit_distance(*d_bak, allow_insert, h_bak, max_k - 1);
			// phase2: evaluate the result
			if(branch_edits <= max_k){
				max_k = branch_edits;
				if(hist) *best_h = *h_bak;
			}
			success |= (branch_edits <= max_k);
		} else eddbgcout << t_arc(path.second) << " is permanent\n";
		delete d_bak;
		if(hist) {
			delete h_bak;
			delete prep_h;
			*hist = *best_h;
			delete best_h;
		}
		if(success) return prep_edits + max_k; else return MY_INFINITY;
	} else { // digraph is transitive
		eddbgcout << "found solution in k = " << hist->size() << " edits ("<<prep_edits<<" preprocessing edits)\n";
		if(hist) delete prep_h;
		return prep_edits;
	}
}

// test whether there is an edit set of size <=k that turns d transitive, return the size of an edit set that is smaller then k
uint test_edit_ex_distance(
		t_dir_graph& d,
		const uint prev_k,
		const bool allow_insert,
		t_edit_ex_history* hist,
		const uint branching_level){
	eddbgcout << "=====> testing digraph (k = "<< prev_k << ") " << d << "\n";
	if(hist) eddbgcout << "edits: " << *hist << "\n";
	if(prev_k >= MY_INFINITY) return MY_INFINITY;
	uint max_k = prev_k;
	uint branch_edits = MY_INFINITY;
	t_edit_ex_history* hist_bak = NULL;
	bool diamond_exists = allow_insert;
	set<t_3path> paths;
	t_vertex v;
	p_diamond diam;

	prep_remove_good(&d);

	if(allow_insert){
		diam = d.get_a_diamond();
		if(diam){
			v = *(set_intersect(d.get_predecessors(diam->second),d.get_successors(diam->first)).begin());
			paths.insert(t_3path(t_arc(diam->first,v),t_arc(v, diam->second)));
			delete diam;
		} else diamond_exists = false;
	}
	/*******************************************************
	 * first of all: check for weakly connected components *
	 *******************************************************/
	set<t_vertex> src_snk;
	t_dir_graph* temp_d = new t_dir_graph(d);
	// if we're calculating transitivity deletion, we can remove all sinks and sources
	// before checking for components. however, we need to insert them back into each
	// component after split up
	if(!diamond_exists){
//		src_snk = set_unite(d.get_sources(), d.get_sinks());
		src_snk = set_substract(d.get_vertices(), d.get_p3_second());
		temp_d->delete_vertices(src_snk);
	}
	t_dir_graph* temp_comp = temp_d->split_component();
	if(temp_d->get_vertices().size()){
		// we have 2 components: temp_d & temp_comp
		if(src_snk.size()){
			// rejoin sources and sinks
			set<t_vertex> v1 = set_unite(temp_d->get_vertices(), src_snk);
			set<t_vertex> v2 = set_unite(temp_comp->get_vertices(), src_snk);
			delete temp_d;
			delete temp_comp;
			temp_d = d.get_induced_subgraph(v1);
			temp_comp = d.get_induced_subgraph(v2);
		}
		prepdbgcout << "found 2 components: " << *temp_d << "\n and " << *temp_comp << "\n";
		uint needed_k = 1;
		if(hist) hist_bak = new t_edit_ex_history(*hist);
		t_dir_graph* temp_d_bak = new t_dir_graph(*temp_d);

		while((needed_k <= max_k) && (test_edit_ex_distance(*temp_d, needed_k, diamond_exists, hist, branching_level+1) >=MY_INFINITY)){
			needed_k++;
			if(hist) *hist = *hist_bak;
			*temp_d = *temp_d_bak;
		}
		delete temp_d_bak;
		delete temp_d;
		if(hist) delete hist_bak;
		// remove from hist and insert back into the graph those arcs between sources and sinks
		max_k -= needed_k;
		prepdbgcout<<"needed "<<needed_k<<" edits for the first component, "<< max_k<< " edits left\n";
		if(hist) prepdbgcout << "history so far: " << *hist << "\n";
		needed_k = test_edit_ex_distance(*temp_comp, max_k, diamond_exists, hist, branching_level+1);
		delete temp_comp;
		return needed_k;
	}
	delete temp_comp; delete temp_d;
	/*******************************
	 * second: apply preprocessing *
	 *******************************/
	if(hist) {prepdbgcout << "hist so far: " << *hist << ", now preprocessing:\n";}
		else {prepdbgcout << "no hist, still preprocessing\n";}
	t_edit_ex_history* prep_hist = NULL;
	if(hist) prep_hist = new t_edit_ex_history();
	// for interleaving the search tree
	if(d.get_vertex_number() > get_kernel_size(max_k))
			max_k = apply_preprocessing(&d, prev_k, diamond_exists, prep_hist);
	uint prep_edits = prev_k - max_k;
	
	prepdbgcout << "after preprocessing: " << d << " with k = " << max_k << ((max_k>=MY_INFINITY)?"(inf)":"");
	if(hist) {prepdbgcout <<" (edits:"<<*hist<< ")\n";} else{ prepdbgcout << " (no hist)\n";}
	if(max_k >= MY_INFINITY) {
		// preprocessing figured that it's not possible with prev_k edits
		if(hist) delete prep_hist;
		return max_k;
	}
	/*******************************************************
	 * third: branch all 3 possibilities to destroy the P3 *
	 *******************************************************/
	paths.clear();
	if(diamond_exists){
		diam = d.get_a_diamond();
		if(diam){
			v = *(set_intersect(d.get_predecessors(diam->second),d.get_successors(diam->first)).begin());
			paths.insert(t_3path(t_arc(diam->first,v),t_arc(v, diam->second)));
			delete diam;
		} else diamond_exists = false;
	}
	p_3path u;
	if(paths.empty()){
		u = d.get_a_p3();
//		u = get_branching_P3(d, d.get_p3(), diamond_exists);
		if(u){
			paths.insert(*u);
			delete u;
		}
		diamond_exists = false;
	}
	if(!paths.empty()) {
		t_3path path = *(paths.begin());
		if(max_k == 0) {
			if(hist) delete prep_hist;
			return MY_INFINITY;
		} // if we dont have any edits left, return infinity
		bool restore1, restore2;
		restore1 = restore2 = false;
		t_dir_graph* d_bak;
		hist_bak = NULL;
		
		// make a copy of the original graph and edit history to work on
		d_bak = new t_dir_graph(d);
		if(hist) {
			hist_bak = new t_edit_ex_history(*hist);
			append_history(hist_bak, prep_hist);
		}
		// === 1. try deleting the first arc ===
		if(d_bak->delete_arc(path.first, true)){
			if(hist_bak) hist_bak->push_back(t_edit_ex(t_edit(ED_DELETE, path.first),""));
			branch_edits = 1 + test_edit_ex_distance(*d_bak, max_k - 1, allow_insert && diamond_exists, hist_bak, branching_level+1);
			if(branch_edits <= max_k) { // found an edit set ! clean up and return the edit set
				if(restore1) d.mark_arc(t_arc(path.first.first, path.second.second), MRK_NONE);
				delete d_bak;
				if(hist){
					*hist = *hist_bak;
					delete hist_bak;
					delete prep_hist;
				}
				return branch_edits + prep_edits; // if we find an edit set that is smaller then max_k, return it
			} else { // an edit set couldnt be found ! mark the arc permanent and restore copies
				d.mark_arc(path.first, MRK_PERMANENT);
				restore2 = true;
				*d_bak = d;
				if(hist){
					*hist_bak = *hist;
					append_history(hist_bak, prep_hist);
				}
			}
		}
		// === 2. try deleting the second arc ===
		if(d_bak->delete_arc(path.second, true)){
			if(hist_bak) hist_bak->push_back(t_edit_ex(t_edit(ED_DELETE, path.second),""));
			branch_edits = 1 + test_edit_ex_distance(*d_bak, max_k - 1, allow_insert && diamond_exists, hist_bak, branching_level+1);
			if(branch_edits <= max_k) { // found an edit set ! clean up and return the edit set
				delete d_bak;
				if(hist){
					*hist = *hist_bak;
					delete hist_bak;
					delete prep_hist;
				}
				return branch_edits + prep_edits; // if we find an edit set that is smaller then max_k, return it
			} else { // an edit set couldnt be found ! mark the arc forbidden and restore copies
				d.mark_arc(path.second, MRK_PERMANENT);
				restore1 = true;
				*d_bak = d;
				if(hist){
					*hist_bak = *hist;
					append_history(hist_bak, prep_hist);
				}
			}
		}
		
		// === 3. insert the missing arc and branch ===
		if(diamond_exists && allow_insert)
			if(d_bak->insert_arc(t_arc(path.first.first, path.second.second),true)){
				if(hist_bak) hist_bak->push_back(t_edit_ex(t_edit(ED_INSERT, t_arc(path.first.first, path.second.second)),""));
				branch_edits = 1 + test_edit_ex_distance(*d_bak, max_k - 1, ED_INSERT, hist_bak, branching_level+1);
			}
		
		if(restore1) d.mark_arc(t_arc(path.first.first, path.second.second), MRK_NONE);
		if(restore2) d.mark_arc(path.first, MRK_NONE);
		
		delete d_bak;

		// found an edit set ! clean up and return the edit set
		if(hist){
			if(branch_edits <= max_k) *hist = *hist_bak;
			delete hist_bak;
			delete prep_hist;
		}
		return (branch_edits <= max_k)?(branch_edits + prep_edits):MY_INFINITY;
	} else { // digraph is transitive ! clean up and return the preprocessing edits
		if(hist){
			append_history(hist, prep_hist);
			delete prep_hist;
		}
		return prep_edits; // number of preprocessing steps
	}
}
// calc the minimum transitive edit distance indirectly, by binary searching k using test_edit_ex_distance
uint get_transitive_min_edit_distance_indirectly(
		t_dir_graph& d,
		const bool allow_insert,
		t_edit_ex_history* best_h){
	uint max_k = MY_INFINITY;
	uint min_k = 1;
	uint k;

	min_k = get_maximal_disjoint_P3(d, get_3paths(d, MY_INFINITY)).size();
	
	t_dir_graph *d2;
	t_edit_ex_history* h = new t_edit_ex_history();
	while(max_k > min_k + 1){
		d2 = new t_dir_graph(d);
		h->clear();
		k = (max_k >= MY_INFINITY)?min_k<<1:(max_k + min_k)>>1;
		eddbgcout << "---------------------------------------------------------------------\n";
		eddbgcout << "checking for k <= " << k << " (k between " << min_k << " and " << max_k << ")\n";
		k = test_edit_ex_distance(*d2, k, allow_insert, h);
		eddbgcout << "got k <= " << k << "\n";
		if(k >= MY_INFINITY) min_k = (max_k >= MY_INFINITY)?min_k<<1:(min_k + max_k)>>1;
			else {
				if(best_h) *best_h = *h;
				max_k = k;
			}
		delete d2;
	}
	delete h;
	return max_k;
}

// calc the minimum transitive edit distance indirectly, by linear search for k using test_edit_ex_distance
uint get_transitive_min_edit_distance_indirectly_linear(
		t_dir_graph& d,
		const bool allow_insert,
		t_edit_ex_history* best_h){
	uint k = get_maximal_disjoint_P3(d, get_3paths(d, MY_INFINITY)).size();
	t_edit_ex_history* h = NULL;
	if(best_h) h = new t_edit_ex_history();
	t_dir_graph* d2 = new t_dir_graph(d);
	while((test_edit_ex_distance(*d2, k, allow_insert, h) >= MY_INFINITY)){
		k++;
		eddbgcout << "---------------------------------------------------------------------\n";
		eddbgcout << "checking for k = " << k << "\n";

		if(h) h->clear();
		*d2 = d;
	}
	if(best_h && h) *best_h = *h;
	delete h;
	delete d2;

	return k;
}

// calc the minimum transitive edit distance by greedy
uint get_transitive_min_edit_distance_greedy(
		t_dir_graph& d,
		const bool allow_insert,
		t_edit_ex_history* best_h,
		map<t_arc, int>* arc_label){
	set<t_vertex> v_set = d.get_vertices();
	set<t_vertex> suc1,suc2,pre1,pre2;
	int modify;

	if(!arc_label){ // initialize arc_label with 0
		arc_label = new map<t_arc, int>();
		for(set<t_vertex>::iterator i = v_set.begin(); i != v_set.end(); i++)
			for(set<t_vertex>::iterator j = v_set.begin(); j != v_set.end(); j++) if(*i != *j){
				// the label of an arc is the gain in chainging this arc
				pre1 = d.get_predecessors(*i); pre1.erase(*j);
				pre2 = d.get_predecessors(*j); pre2.erase(*i);
				suc1 = d.get_successors(*i); suc1.erase(*j);
				suc2 = d.get_successors(*j); suc2.erase(*i);
				modify = set_substract(pre1, pre2).size()
						+set_substract(suc2, suc1).size()
						-set_intersect(suc1, pre2).size();
				if(!d.contains_arc(t_arc(*i,*j))) modify = -modify;
//				greeddbgcout << "marking (" << *i << "," << *j << ") with "<< modify << " [pre(" << *i << ")=" << pre1 << " suc(" << *i << ")=" << suc1 << " pre(" << *j << ")=" << pre2 << " suc(" << *j << ")=" << suc2 << "\n";
				arc_label->insert(pair<t_arc, int>(t_arc(*i,*j),modify));
			}
	}

	greeddbgcout << "arc_label:\n ";
	for(map<t_arc, int>::const_iterator i = arc_label->begin(); i != arc_label->end(); i++)
		if(i->second > 0) greeddbgcout << "(" << i->first << "," << i->second << "),";
	greeddbgcout << "\n";

	map<t_arc, int>::iterator max_arc = arc_label->end();
	for(map<t_arc, int>::iterator i = arc_label->begin(); i != arc_label->end(); i++)
//		if(!d.is_permanent(i->first)){
//			if(!d.is_forbidden(i->first)){
				if((d.contains_arc(i->first)) || (i->second > 1))
					if((max_arc == arc_label->end()) || (i->second > max_arc->second)) max_arc = i;
//			} else greeddbgcout << i->first << " is forbidden\n";
//		} else greeddbgcout << i->first << "is permanent\n";
	if(max_arc == arc_label->end()){
		delete arc_label;
		return 0;
	}
	if(max_arc != arc_label->end()) greeddbgcout << "max_arc: " << max_arc->first << " (" << max_arc->second << ")\n";
	if(max_arc->second > 0){
		t_vertex v;
		t_edit_ex e;
		t_edit_ex_history::iterator same_arc;
		bool opposite_in_h = false;
		stringstream explain;
		short exist_1, exist_2;
		t_edit_ex_history* h;
		if(best_h) h = best_h; else h = new t_edit_ex_history();

		same_arc = h->end();
		for(t_edit_ex_history::iterator i = h->begin(); i != h->end(); i++)
			if((i->first.second == max_arc->first)){
				same_arc = i;
				break;
			}

		if((same_arc == h->end()) || ((same_arc->first.first == ED_DELETE) == (d.contains_arc(max_arc->first)))){
			opposite_in_h = false;
			explain.str("");
			explain << "G(" << max_arc->second << ")";
			if(d.contains_arc(max_arc->first)) {
				d.delete_arc(max_arc->first, false);	
				h->push_back(t_edit_ex(t_edit(ED_DELETE, max_arc->first),explain.str()));	
				modify = 1;
				greeddbgcout << max_arc->first << " is deleted\n";
				
			} else {
				d.insert_arc(max_arc->first, false);
				h->push_back(t_edit_ex(t_edit(ED_INSERT, max_arc->first),explain.str()));
				modify = -1;
				greeddbgcout << max_arc->first << " is inserted\n";
			}
			(*arc_label)[max_arc->first] = -(*arc_label)[max_arc->first];
			// update P3 destruction ratings
			for(set<t_vertex>::const_iterator u = v_set.begin(); u != v_set.end(); u++)
				if((*u != max_arc->first.first) && (*u != max_arc->first.second)){
					if(d.contains_arc(t_arc(*u, max_arc->first.first))) exist_1 = 1; else exist_1 = -1;
					if(d.contains_arc(t_arc(*u, max_arc->first.second))) exist_2 = 1; else exist_2 = -1;
						(*arc_label)[t_arc(*u,max_arc->first.first)] += ((exist_2 - 1)>>1) * exist_1 * modify;
						(*arc_label)[t_arc(*u,max_arc->first.second)] += ((exist_1 + 1)>>1) * exist_2 * modify;
					if(d.contains_arc(t_arc(max_arc->first.first, *u))) exist_1 = 1; else exist_1 = -1;
					if(d.contains_arc(t_arc(max_arc->first.second, *u))) exist_2 = 1; else exist_2 = -1;
						(*arc_label)[t_arc(max_arc->first.first, *u)] += ((exist_2 + 1)>>1) * exist_1 * modify;
						(*arc_label)[t_arc(max_arc->first.second, *u)] += ((exist_1 - 1)>>1) * exist_2 * modify;
					if(d.contains_arc(t_arc(max_arc->first.first, *u))) exist_1 = 1; else exist_1 = -1;
					if(d.contains_arc(t_arc(*u, max_arc->first.second))) exist_2 = 1; else exist_2 = -1;
						(*arc_label)[t_arc(max_arc->first.first, *u)] += ((exist_2 + 1)>>1) * exist_1 * modify;
						(*arc_label)[t_arc(*u, max_arc->first.second)] += ((exist_1 + 1)>>1) * exist_2 * modify;
			}
		} else {
			cout << "found opposite arc: " << *same_arc << "\n";
			// delete same_arc from history
			h->erase(same_arc);
			opposite_in_h = true;
		}
		uint res = get_transitive_min_edit_distance_greedy(d, allow_insert, h, arc_label);
		if(!best_h) delete h;
		if(opposite_in_h) return res - 1; else return res + 1;
	} else {
		delete arc_label;
		return 0;
	}
}


void transitive_closure(t_dir_graph* d, t_edit_ex_history* hist){
	p_3path path;
	while((path = d->get_a_p3())){
		if(d->insert_arc(t_arc(path->first.first, path->second.second)))
			if(hist) 
				hist->push_back(t_edit_ex(t_edit(ED_INSERT,
							t_arc(path->first.first, path->second.second)),""));
	}
}




#endif

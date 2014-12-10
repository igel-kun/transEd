/***************************************************
 * editing.h
  * by M.Weller
 **************************************************/
 
#ifndef editing_h
#define editing_h

#include "digraph.h"
#include <vector>

#define ED_DEBUG 0
#define PREP_DEBUG 0
#define GREED_DEBUG 0
#define eddbgcout if(ED_DEBUG) cout
#define prepdbgcout if(PREP_DEBUG) cout
#define greeddbgcout if(GREED_DEBUG) cout

#define ED_INSERT true
#define ED_DELETE false

#define INTERLEAVING_CONSTANT 0.01
#define SPLIT_SEPERATOR (char)'\''
using namespace std;

typedef pair<bool, t_arc> t_edit; // true = insert
typedef pair<t_edit, string> t_edit_ex; // string for explanation
typedef vector<t_edit> t_edit_history;
typedef vector<t_edit_ex> t_edit_ex_history;


ostream& operator<<(ostream& os, const t_3path& p);
ostream& operator<<(ostream& os, const set<t_3path>& sp);
ostream& operator<<(ostream& os, const t_edit& e);
ostream& operator<<(ostream& os, const t_edit_ex& e_ex);
ostream& operator<<(ostream& os, const t_edit_ex_history& h);

inline uint get_kernel_size(const uint k){
	return (uint)(INTERLEAVING_CONSTANT * (double)k * (double)(k+2));
}

void append_history(t_edit_ex_history* h1, const t_edit_ex_history* h2);
set<t_3path> get_maximal_disjoint_P3(const t_dir_graph& d, const set<t_3path>& paths, const bool allow_insert = true);

// apply preprocessing rule2: remove doublearc-parts where possible
void prep_remove_double(t_dir_graph* d);
// apply preprocessing rule: remove all good vertices
void prep_remove_good(t_dir_graph* d);
// apply preprocessing rule: process high-degree rules
uint prep_high_degree(
		t_dir_graph* d,
		const uint prev_k,
		const bool allow_insert = ED_INSERT,
		t_edit_ex_history* h = NULL);
// apply preprocessing rule: cut the endings and return the new k
uint prep_cut_endings(
		t_dir_graph* d,
		const uint prev_k,
		const bool allow_insert = ED_INSERT,
		t_edit_ex_history* h = NULL);
//  apply all preprocessing rules
uint apply_preprocessing(
		t_dir_graph* d,
		const uint prev_k,
		const bool allow_insert = ED_INSERT,
		t_edit_ex_history* h = NULL);

// check if the given triple is a P3
bool is_3path(const t_dir_graph& d,const t_vertex& u,const t_vertex& v,const t_vertex& w);
// return a set of induced P3 that is in d containing at max max_p3 P3
set<t_3path> get_3paths(const t_dir_graph& d, const uint max_p3 = 1);
// returns a set of (>min_n-1)-diamonds represented as one of its 3pathes in d of size at most max_diamonds
set<t_3path> get_diamond(const t_dir_graph& d, uint max_diamonds = 1, uint min_n = 2);

// calc the minimum transitive edit distance  directly
uint get_transitive_min_edit_distance(
		t_dir_graph& d, 
		const bool allow_insert = ED_INSERT,
		t_edit_ex_history* hist = NULL,
		const uint prev_k = MY_INFINITY-2);

// test whether there is an edit set of size <=k that turns d transitive, return the size of an edit set that is smaller then k
uint test_edit_ex_distance(
		t_dir_graph& d,
		const uint prev_k,
		const bool allow_insert = ED_INSERT,
		t_edit_ex_history* best_h = NULL,
		const uint branching_level = 0);

// calc the minimum transitive edit distance indirectly, by binary searching k using test_edit_ex_distance
uint get_transitive_min_edit_distance_indirectly(
		t_dir_graph& d,
		const bool allow_insert = ED_INSERT,
		t_edit_ex_history* best_h = NULL);

// calc the minimum transitive edit distance indirectly, by binary searching k using test_edit_ex_distance
uint get_transitive_min_edit_distance_indirectly_linear(
		t_dir_graph& d,
		const bool allow_insert = ED_INSERT,
		t_edit_ex_history* best_h = NULL);
// calc the minimum transitive edit distance by greedy
uint get_transitive_min_edit_distance_greedy(
		t_dir_graph& d,
		const bool allow_insert = ED_INSERT,
		t_edit_ex_history* best_h = NULL,
		map<t_arc, int>* arc_label = NULL);

// calc the transitive closure of a digraph
void transitive_closure(t_dir_graph* d, t_edit_ex_history* hist = NULL);



#endif

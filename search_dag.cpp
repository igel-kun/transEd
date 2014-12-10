#include <iostream>
#include "digraph.h"
#include "editing.h"


void build_skeleton(t_dir_grah* D, set<t_vertex>* available){
	D->insert_vertex("a1");
	D->insert_vertex("a2");
	D->insert_vertex("a3");
	D->insert_vertex("a4");
	available->insert("a2");
	available->insert("a3");
	
	D->insert_vertex("b1");
	D->insert_vertex("b2");
	D->insert_vertex("b3");
	D->insert_vertex("b4");	
	available->insert("b2");
	available->insert("b3");
	
	D->insert_vertex("c1");
	D->insert_vertex("c2");
	D->insert_vertex("c3");
	D->insert_vertex("c4");
	available->insert("c2");
	available->insert("c3");
}

void set_var(t_dir_graph* D, const string var, const bool value){
	t_arc a1(var+"1", var+"2");
	t_arc a2(var+"2", var+"3");
	t_arc a3(var+"3", var+"4");
	D->mark_arc(a1,MRK_NONE);D->mark_arc(a2,MRK_NONE);D->mark_arc(a3,MRK_NONE);
	if(value){
		D->insert_arc(a1,true);
		D->delete_arc(a2,true);
		D->insert_arc(a3,true);
	} else {
		D->delete_arc(a1,true);
		D->insert_arc(a2,true);
		D->delete_arc(a3,true);
	}
}

bool test_graph(const t_dir_graph* D, const uint a){
}

int main(int argc, char** argv){
	uint k;
	t_dir_graph* D = new t_dir_graph();
	set<t_vertex> available;
	
	ostringstream vconstruct;
	t_vertex new_v;
	uint count = 0;
	uint poss;
	
	build_skeleton(D, &available);
	while(1){
		poss = 1 << (available.size() * availalbe.size());
		t_dir_graph *Dprime = new t_dir_graph(*D);
		for(uint j = 0; j < poss; j++){
			// build the arrangement and test it
			for(set<t_vertex>::const_iterator u = available.begin(); u != available.end(); u++){
				for(set<t_vertex>::const_iterator v = available.begin(); v != available.end(); v++){
				
			
			// test phase
			for(uint a = 0; a < 8 ; a++) test_graph(D,a);
		}
		// if no arrangement could be found, take an additional vertex
		vconstruct << "x" << (count++);
		new_v = vconstruct.str();
		available.insert(new_v);
		D->insert_vertex(new_v);
	}

	t_edit_history* best_h = new t_edit_history();
	the_time = time(NULL);
	k = get_transitive_min_edit_distance_indirectly(*Dprime, best_h);
	cout << k << " edits neccessary: " << *best_h << " (time used: " << the_time - time(NULL) <<")\n";
	
	cout << "\nnow directly:\n";
	delete best_h;best_h = new t_edit_history();
	the_time = time(NULL);
	k = get_transitive_min_edit_distance(*D, best_h);
	cout << k << " edits neccessary: " << *best_h << " (time used: " << the_time - time(NULL) <<")\n";
	delete best_h;
	delete D;
}

#include <ctime>
#include <iostream>
#include "digraph.h"
#include "editing.h"

#include <math.h>
void cap_k(t_dir_graph* d, const uint max_k){
	t_arc a;
	set<t_vertex> vertices = d->get_vertices();
	uint r;
	set<t_vertex>::iterator j;
	t_vertex v;

	transitive_closure(d);
	srand(time(NULL));
	for(uint i = 0; i < max_k; i++){
		do{
			// get a random arc
			r = (uint)floor((((rand()+1)*(double)vertices.size()))/(double)RAND_MAX);
			j = vertices.begin();
			while(r){
				r--;
				j++;
			}
			v = *j;

			r = (uint)floor((((rand()+1)*(double)vertices.size()))/(double)RAND_MAX);
			j = vertices.begin();
			while(r){
				r--;
				j++;
			}
		} while(v == *j);
		a = t_arc(v,*j);
		if(d->contains_arc(a)) d->delete_arc(a); else d->insert_arc(a);
	}
}


int main(int argc, char** argv){
	uint k;
	t_dir_graph* D;
	t_edit_ex_history* best_h;
	uint the_time;

	if(argc > 3){
		cout << "generating digraph ("<<argv[1]<<" vertices, edge probability "<<argv[2]<<"%, max k="<< argv[3]<<")\n";
		D = generate_random_digraph(atoi(argv[1]), (double)atoi(argv[2])/(double)100);
		cout << "capping k...\n";
		cap_k(D,atoi(argv[3]));
	} else if(argc > 2) {
		cout << "generating digraph ("<<argv[1]<<" vertices, edge probability "<<argv[2]<<"%)\n";
		D = generate_random_digraph(atoi(argv[1]), (double)atoi(argv[2])/(double)100);
	} else if(argc > 1){
		char* filename = argv[1];
		printf("reading file %s\n", filename);
		D = new t_dir_graph();
		D->read_from_file(filename);
	} else {
		//printf("input a digraph:\n");
		//D = new t_dir_graph();
		//D->read_from_stream(cin);
    cout << "syntax: " << argv[0] << " (file| <#vertices> <edge probability> [maximum k])" <<endl;
    exit(1);
	}
	cout << "clocks per second: " << CLOCKS_PER_SEC << "\n";
//	cout << "removing good vertices...\n";
//	prep_remove_good(D);
	cout << "graph: " << *D << "\n";

	t_dir_graph *Dprime = new t_dir_graph(*D);
	cout << "number of P3: " << D->get_p3().size() << " ("<< get_maximal_disjoint_P3(*D, D->get_p3(), true).size() << " disjoint)\n";
	
	cout << "\ncomputing greedy heuristic:\n";
	best_h = new t_edit_ex_history();
	the_time = clock();
	k = get_transitive_min_edit_distance_greedy(*Dprime, ED_INSERT, best_h);
	the_time = clock() - the_time;
	cout << k << " edits neccessary: " << *best_h << " (time used: " << ((double)the_time)/CLOCKS_PER_SEC <<")\n";
	delete best_h;
	delete Dprime;

	Dprime = new t_dir_graph(*D);
	cout << "\ncomputing exact solution:\n";
	best_h = new t_edit_ex_history();
	the_time = clock();
	k = get_transitive_min_edit_distance_indirectly_linear(*Dprime, ED_INSERT, best_h);
	the_time = clock() - the_time;
	cout << k << " edits neccessary: " << *best_h << " (time used: " << ((double)the_time)/CLOCKS_PER_SEC <<")\n";
	delete best_h;
	delete Dprime;
/*
	Dprime = new t_dir_graph(*D);
	cout << "\nnow directly:\n";
	best_h = new t_edit_ex_history();
	the_time = clock();
	k = get_transitive_min_edit_distance(*Dprime, ED_INSERT, best_h, k); // greedy is upper bound
	the_time = clock() - the_time;
	cout << k << " edits neccessary: " << *best_h << " (time used: " << ((double)the_time)/CLOCKS_PER_SEC <<")\n";
	delete best_h;
	delete Dprime;
*/
	delete D;
}

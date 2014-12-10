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
	for(uint i = 0; i < max_k; i++){
		do{
			// get a random arc
			r = (uint)floor(((rand()+1)*vertices.size())/RAND_MAX);
			j = vertices.begin();
			while(r){
				r--;
				j++;
			}
			v = *j;

			r = (uint)floor(((rand()+1)*vertices.size())/RAND_MAX);
			j = vertices.begin();
			while(r){
				r--;
				j++;
			}
			a = t_arc(v,*j);
		} while(v == *j);

		if(d->contains_arc(a)) d->delete_arc(a); else d->insert_arc(a);
	}
}


t_dir_graph* get_graph_by_number(const uint nr){
	t_dir_graph *d = new t_dir_graph();
	const int n = 6; // remember: we need n*(n-1) bits in nr
	char i[2];
	char j[2];

	i[1]=0;
	j[1]=0;
	for(i[0] = 48; i[0]-48 < n; i[0]++) d->insert_vertex(t_vertex(i));
	for(i[0] = 48; i[0]-48 < n; i[0]++)
		for(j[0] = 48; j[0]-48 < n; j[0]++) if(i[0]!=j[0])
			if((nr >> (n*(i[0]-48)+j[0]-48-1)) & 1)
				d->insert_arc(t_arc(t_vertex(i),t_vertex(j)));
	return d;
}


int main(int argc, char** argv){
	uint k;
	t_dir_graph* D;
	t_dir_graph *Dprime;
	uint dir_the_time;
	uint dir_k;
	p_3path p3;
	uint max_n[36];
	srand(time(NULL));

	for(uint i = 0; i < 36; i++) max_n[i]=0;

//	for(uint horst = 0; horst < 100000; horst++)
	while(true)
	for(uint i = 1; i < 100; i++){
//	for(uint i = 0; i < 0xffffffff; i++){
		D = generate_random_digraph(18,((double)i)/100);
//		D = get_graph_by_number(i);
//		cap_k(D,9);
		p3 = get_3path(*D);
		if(p3){
			delete p3;
			Dprime = new t_dir_graph(*D);
//			cout << "graph: " << *D << "\n";
	
			dir_the_time = time(NULL);
//			dir_k = get_transitive_min_edit_distance_indirectly_linear(*D);
			dir_k = test_edit_ex_distance(*D, 8);
			dir_the_time = time(NULL) - dir_the_time;
//			cout << dir_k << " edits neccessary: " << *dir_best_h << " (time used: " << dir_the_time <<")\n";
			if(dir_k < MY_INFINITY){
				delete D;
				D = new t_dir_graph(*Dprime);
				prep_remove_good(D);

				k = apply_preprocessing(Dprime, dir_k);
				if((k == dir_k) && (D->get_vertices().size() > max_n[k])) {
					max_n[k]=D->get_vertices().size();
					cout << "new max_n["<<k<<"]="<<max_n[k]<<": " << *D << "\n";
				}
			}
			delete Dprime;
		}
		delete D;
	}
}

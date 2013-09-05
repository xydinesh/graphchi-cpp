/**
 * @file
 * @author  Aapo Kyrola <akyrola@cs.cmu.edu>
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright [2012] [Aapo Kyrola, Guy Blelloch, Carlos Guestrin / Carnegie Mellon University]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 
 *
 * @section DESCRIPTION
 *
 * Template for GraphChi applications. To create a new application, duplicate
 * this template.
 */



#include <string>
#include <vector>
#include "graphchi_basic_includes.hpp"
#include "engine/dynamic_graphs/graphchi_dynamicgraph_engine.hpp"
#include "engine/auxdata/degree_data.hpp"
#include "preprocessing/util/orderbydegree.hpp"

#define NBRSIZE 500
using namespace graphchi;

int high_degree_node = 20;
int grabbed_edges = 0;
FILE *out_file;
long long count = 1;

struct dense_adj {
  int count;
  vid_t * adjlist;

  dense_adj() { adjlist = NULL; }
  dense_adj(int _count, vid_t * _adjlist) : count(_count), adjlist(_adjlist) {
  }

};


// This is used for keeping in-memory
class adjlist_container {
  std::vector<dense_adj> adjs;
  //mutex m;
public:
  vid_t pivot_st, pivot_en;
  std::vector<bool> high_degree_nodes;

  adjlist_container() {
    pivot_st = 0; //start pivor on item nodes (excluding user nodes)
    pivot_en = 0;
  }

  void clear() {
    for(std::vector<dense_adj>::iterator it=adjs.begin(); it != adjs.end(); ++it) {
      if (it->adjlist != NULL) {
        free(it->adjlist);
        it->adjlist = NULL;
      }
    }
    adjs.clear();
    high_degree_nodes.clear();
    pivot_st = pivot_en;
  }

  /** 
   * Extend the interval of pivot vertices to en.
   */
  void extend_pivotrange(vid_t en) {
    assert(en>=pivot_en);
    pivot_en = en; 
    adjs.resize(pivot_en - pivot_st);
    high_degree_nodes.resize(pivot_en - pivot_st);
  }

  int load_edges_into_memory(graphchi_vertex<uint32_t, uint32_t> &v) {

    int num_edges = v.num_edges();
    //not enough user rated this item, we don't need to compare to it
    if (num_edges > high_degree_node){
      high_degree_nodes[v.id()] = true;
      return 0;
    }
       
    high_degree_nodes[v.id()] = false;
    // Count how many neighbors have larger id than v
    dense_adj dadj = dense_adj(num_edges, (vid_t*) calloc(sizeof(vid_t), num_edges));
    for(int i=0; i<num_edges; i++) {
      dadj.adjlist[i] = v.edge(i)->vertex_id();
    }

    //std::sort(dadj.adjlist, dadj.adjlist + num_edges);
    adjs[v.id() - pivot_st] = dadj;
    assert(v.id() - pivot_st < adjs.size());
    __sync_add_and_fetch(&grabbed_edges, num_edges /*edges_to_larger_id*/);
    return num_edges;
  }

  void print_fof(vid_t v) {
    dense_adj dv = adjs[v];
    //int rc = fprintf(out_file, "%d:%d \n", v, dv.count);
    if (high_degree_nodes[v] == true) {
      // rc = fprintf(out_file, "high degree node > %d \n", high_degree_node);
      return;
    }
    
    for (int i = 0; i < dv.count; i++) {
      // rc = fprintf(out_file, "\t %d| ", dv.adjlist[i]); 
      // if (rc <= 0) {
      // 	perror("Failed to write output\n");
      // 	logstream(LOG_FATAL)<<"Failed to write output to: .out" << std::endl;
      // }
      
      if (high_degree_nodes[dv.adjlist[i]] == true) {
	// rc = fprintf(out_file, "high degree node > %d \n", high_degree_node);
	continue;
      }
      
      dense_adj ddv = adjs[dv.adjlist[i]];
      for (int j = 0; j < ddv.count; j++) {
	// rc = fprintf(out_file, "%d, ", ddv.adjlist[j]); 
	// if (rc <= 0) {
	//   perror("Failed to write output\n");
	//   logstream(LOG_FATAL)<<"Failed to write output to: .out" << std::endl;
	// }
      }

      // rc = fprintf(out_file, "\n");
      // if (rc <= 0) {
      // 	perror("Failed to write output\n");
      // 	logstream(LOG_FATAL)<<"Failed to write output to: .out" << std::endl;
      // }
    }

    // rc = fprintf(out_file, "\n");
    // if (rc <= 0) {
    //   perror("Failed to write output\n");
    //   logstream(LOG_FATAL)<<"Failed to write output to: .out" << std::endl;
    // }
  }

  inline bool is_pivot(vid_t vid) {
    return (vid >= pivot_st) && (vid < pivot_en);
  }
};


adjlist_container * adjcontainer;
/**
 * Type definitions. Remember to create suitable graph shards using the
 * Sharder-program. 
 */
typedef unsigned int VertexDataType;
typedef unsigned int EdgeDataType;


/*
 * Class for writing the output number of triangles for each node
 */
class OutputVertexCallback : public VCallback<VertexDataType> {
public:
  virtual void callback(vid_t vertex_id, VertexDataType &value) {
    if (value > 0)
      std::cout << vertex_id << " " << value << std::endl;
  }
};


/**
 * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
 * class. The main logic is usually in the update function.
 */
struct DegreeDist : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
 
  /**
   *  Vertex update function.
   */
  void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {

    if ((gcontext.iteration % 2) == 0) {
      /* On first iteration, initialize vertex (and its edges). This is usually required, because
         on each run, GraphChi will modify the data files. To start from scratch, it is easiest
         do initialize the program in code. Alternatively, you can keep a copy of initial data files. */ 
      //vertex.set_data(vertex.id());
      if (count % 10000000 == 0)
	printf("even vertex : %lld\n", count);
      adjcontainer->load_edges_into_memory(vertex);
    } else {
      if (count % 10000000 == 0)
	printf("odd vertex id: %lld\n", count);
      adjcontainer->print_fof(vertex.id());
    }
    __sync_add_and_fetch(&count, 1);
  }
    
  /**
   * Called before an iteration starts.
   */
  void before_iteration(int iteration, graphchi_context &gcontext) {
    printf("working on iteration %d\n", iteration);
    logstream(LOG_DEBUG) << "pivot_st: " << adjcontainer->pivot_st << " pivot_en: " << adjcontainer->pivot_en << std::endl;;
    gcontext.scheduler->remove_tasks(0, gcontext.nvertices - 1);
    if (gcontext.iteration % 2 == 0) {
      // Schedule vertices that were pivots on last iteration, so they can
      // keep count of the triangles counted by their lower id neighbros.
      for(vid_t i=adjcontainer->pivot_st; i < adjcontainer->pivot_en; i++) {
	gcontext.scheduler->add_task(i); 
      }
      grabbed_edges = 0;
      adjcontainer->clear();
    } else {
      // Schedule everything that has id < pivot
      logstream(LOG_INFO) << "Now pivots: " << adjcontainer->pivot_st << " " << adjcontainer->pivot_en << std::endl;
      for(vid_t i=0; i < gcontext.nvertices; i++) {
	if (i < adjcontainer->pivot_en) { 
	  gcontext.scheduler->add_task(i); 
	}
      }
    }
  }
    
  /**
   * Called after an iteration has finished.
   */
  void after_iteration(int iteration, graphchi_context &gcontext) {
  }
    
  /**
   * Called before an execution interval is started.
   */
  void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) { 
    printf("starting exec interval %d:%d\n", window_st, window_en);
    if (gcontext.iteration % 2 == 0) {
      if (adjcontainer->pivot_st <= window_en) {
	//size_t max_grab_edges = get_option_long("membudget_mb", 1024) * 1024 * 1024 / 8;
	//if (grabbed_edges < max_grab_edges * 0.8) {
	logstream(LOG_DEBUG) << "Window init, grabbed: " << grabbed_edges << " edges" << std::endl;
	for(vid_t vid=window_st; vid <= window_en; vid++) {
	  gcontext.scheduler->add_task(vid);
	}

	adjcontainer->extend_pivotrange(window_en + 1);
	if (window_en == gcontext.nvertices) {
	  // Last iteration needed for collecting last triangle counts
	  gcontext.set_last_iteration(gcontext.iteration + 2);                    
	}
	//} else {
	//std::cout << "Too many edges, already grabbed: " << grabbed_edges << std::endl;
	//}
      }
    }
  }
    
  /**
   * Called after an execution interval has finished.
   */
  void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
  }
    
};

int main(int argc, const char ** argv) {
  /*
   * Obtained from google groups
   * execthreads = the number of threads used to execute the update function in parallels
   * These are harder to exactly describe: niothreads = num of threads used for asynchronous I/O. 
   * In my experience 2-4 is best. loadthreads = num of threads used in parallel loading of the shards. 2-4 is again good. 
   */ 

  /* GraphChi initialization will read the command line 
     arguments and the configuration file. */
  graphchi_init(argc, argv);
    
  /* Metrics object for keeping track of performance counters
     and other information. Currently required. */
  metrics m("degree-dist");
    
  /* Basic arguments for application */
  std::string filename = get_option_string("file");  // Base filename
  int niters           = get_option_int("niters", 2); // Number of iterations
  //    bool scheduler       = get_option_int("scheduler", 0); //
  //    Whether to use selective scheduling
  high_degree_node           = get_option_int("highdegree", 20); // Number of iterations
  bool scheduler = true;
    
  /* Detect the number of shards or preprocess an input to create them */
  int nshards = convert_if_notexists<EdgeDataType>(filename, 
						   get_option_string("nshards", "auto"));

  //order_by_degree<EdgeDataType>(filename, nshards, m);
  //initialize data structure which saves a subset of the items (pivots) in memory
  adjcontainer = new adjlist_container();
  //array for marking which items are conected to the pivot items via users.
  //high_degree_nodes = new bool[N];
  
  /* Run */
  DegreeDist program;
  graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 
  //engine.set_enable_deterministic_parallelism(false);
  // Low memory budget is required to prevent swapping as triangle counting
  // uses more memory than standard GraphChi apps.
  engine.set_membudget_mb(std::min(get_option_int("membudget_mb", 1024), 1024)); 


  // out_file = fopen((const char*)"/backups/data/fof.txt", "w");
  // if (out_file == NULL) {
  //   perror("fopen failed");
  //   logstream(LOG_FATAL) <<" Failed to open file" << "fof.txt" << std::endl;
  // }

  engine.run(program, niters);
    
  /* Report execution metrics */
  metrics_report(m);

  //OutputVertexCallback callback;
  // foreach_vertices<VertexDataType>(filename + "_degord", 0, engine.num_vertices(), callback);
  return 0;
}

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

#define NBRSIZE 500
using namespace graphchi;

/**
 * Unlike in connected components, we need
 * to ensure that neighbors do not overwrite each
 * others values. This is achieved by keeping two values
 * in an edge. In this struct, smaller_one is the id of the
 * vertex that has smaller id, and larger_one the others.
 * This complexity is due to us ignoring the direction of an edge.
 */
struct edge_label {
  std::vector<vid_t> smaller_one;
  std::vector<vid_t> larger_one;
};

std::vector<vid_t> & neighbor_label(edge_label & bidir, vid_t myid, vid_t nbid) {
  if (myid < nbid) {
    return bidir.larger_one;
  } else {
    return bidir.smaller_one;
  }
}

std::vector<vid_t> & my_label(edge_label & bidir, vid_t myid, vid_t nbid) {
  if (myid < nbid) {
    return bidir.smaller_one;
  } else {
    return bidir.larger_one;
  }
}

/**
 * Type definitions. Remember to create suitable graph shards using the
 * Sharder-program. 
 */
typedef vid_t VertexDataType;
typedef edge_label EdgeDataType;

/**
 * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
 * class. The main logic is usually in the update function.
 */
struct DegreeDist : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
 
  /**
   *  Vertex update function.
   */
  void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {

    if (gcontext.iteration == 0) {
      /* On first iteration, initialize vertex (and its edges). This is usually required, because
         on each run, GraphChi will modify the data files. To start from scratch, it is easiest
         do initialize the program in code. Alternatively, you can keep a copy of initial data files. */
      //vertex.set_data(vertex.id());
      //gcontext.scheduler->add_task(vertex.id()); 
      std::vector<vid_t> nbrs;
      for(int i=0; i < vertex.num_edges(); i++) {
        nbrs.push_back(vertex.edge(i)->vertex_id());
      }
      
      for(int i=0; i < vertex.num_edges(); i++) {
        edge_label elabel = vertex.edge(i)->get_data();
        my_label(elabel, vertex.id(), vertex.edge(i)->vertex_id()) = nbrs;
        vertex.edge(i)->set_data(elabel);
      }
      //vertex.set_data(vertex.id());
    } else {
      /* Do computation */ 
      /* Loop over in-edges (example) */
      // for(int i=0; i < vertex.num_inedges(); i++) {
      //     // Do something
      // //    value += vertex.inedge(i).get_data();
      // }
            
      /* Loop over out-edges (example) */
      // for(int i=0; i < vertex.num_outedges(); i++) {
      //     // Do something
      //     // vertex.outedge(i).set_data(x)
      // }
            
      /* Loop over all edges (ignore direction) */
      printf("%d|%d : \n", vertex.id(), vertex.num_edges());
      for(int i=0; i < vertex.num_edges(); i++) {
        edge_label elabel = vertex.edge(i)->get_data();
        std::vector<vid_t> fof = neighbor_label(elabel, vertex.id(), vertex.edge(i)->vertex_id());
        int len = fof.size();
        printf("\t %d|%d\n", vertex.edge(i)->vertex_id(), len);
        for (int i = 0; i < len; i++) {
          printf("%d|", fof[i]);
        }
        printf("\n");
      }
      printf("\n");

      // for(int i=0; i < vertex.num_edges(); i++) {
      //   edge_label edgelabel = vertex.edge(i)->get_data();
      //   my_label(edgelabel, vertex.id(), vertex.edge(i)->vertex_id()) = friendsvector;
      // }


                      
      //v.set_data(new_value);
    }
  }
    
  /**
   * Called before an iteration starts.
   */
  void before_iteration(int iteration, graphchi_context &gcontext) {
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
  }
    
  /**
   * Called after an execution interval has finished.
   */
  void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
  }
    
};

int main(int argc, const char ** argv) {
  /* GraphChi initialization will read the command line 
     arguments and the configuration file. */
  graphchi_init(argc, argv);
    
  /* Metrics object for keeping track of performance counters
     and other information. Currently required. */
  metrics m("degree-dist");
    
  /* Basic arguments for application */
  std::string filename = get_option_string("file");  // Base filename
  int niters           = get_option_int("niters", 4); // Number of iterations
  //    bool scheduler       = get_option_int("scheduler", 0); //
  //    Whether to use selective scheduling
  bool scheduler = false;
    
  /* Detect the number of shards or preprocess an input to create them */
  int nshards          = convert_if_notexists<EdgeDataType>(filename, 
                                                            get_option_string("nshards", "auto"));
    
  /* Run */
  DegreeDist program;
  graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 
  engine.run(program, niters);
    
  /* Report execution metrics */
  metrics_report(m);
  return 0;
} 

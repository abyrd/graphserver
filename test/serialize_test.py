#! python
from graphserver.core import Graph, State, Street
from graphserver.graphdb import GraphDatabase
import random
import time

def compare_contents(g1, g2) :
    global n_real_errors
    global n_discrepancies
    
    vertices = [v for v in g1.vertices]
    vlabels  = [v.label for v in g1.vertices]
    vlabels.sort()
    
    print "--- Checking that vertex label lists are identical...",
    print "(%d elements)" % len(vlabels)
    vlabels2 = [v.label for v in g2.vertices]
    vlabels2.sort()
    if vlabels != vlabels2 :
        print "Vertices are different - skipping tests for this graph.\n"
        n_real_errors += 1
        n_discrepancies += 1
        return
    
    result = True
    print "--- Checking that vertex payloads are identical..."
    print "(%d elements)" % len(vertices)
    for v1 in vertices :
        v2 = g2.get_vertex(v1.label)
        if v1.payload :
           vp1 = v1.payload.to_xml() 
        else :
           vp1 = "NONE"
        if v2.payload :
           vp2 = v2.payload.to_xml()
        else :
           vp2 = "NONE"
        if not vp1 == vp2 :
            print "discrepancy at :", v1.label
            n_discrepancies += 1
            vw1 = v1.payload.weight
            vt1 = v1.payload.time
            vw2 = v2.payload.weight
            vt2 = v2.payload.time
            if vw1 == vw2 and vt1 == vt2 :
                print "but times and weights are identical."
            elif vt1 == vt2 :
                print "but times are identical."
            else :
                result = False
                print vp1
                print "versus"
                print vp2
                print
                n_real_errors += 1
                            
    print "--- Checking that edge lists and edges are identical..."
    n_vertices = len(vertices)
    n = 0
    for v1 in vertices :
        n += 1
        if n % 10000 == 0 : print "%d / %d" % (n, n_vertices)
        v2 = g2.get_vertex(v1.label)
        ep1 = ["%s to %s %s" % (e.from_v.label, e.to_v.label, e.payload.to_xml()) for e in v1.outgoing] 
        ep2 = ["%s to %s %s" % (e.from_v.label, e.to_v.label, e.payload.to_xml()) for e in v2.outgoing] 
        ep1.sort() 
        ep2.sort() 
        if ep1 != ep2 :
            print "discrepancy at:", v1.label
            n_discrepancies += 1
            if not result : # i.e. if there were differences in the payload times
                print ep1
                print "versus"
                print ep2
                print
                n_real_errors += 1
    
    print
            
# >>> time.ctime(1253700000)
# 'Wed Sep 23 12:00:00 2009'

# >>> time.ctime(1253760000)
# 'Thu Sep 24 04:40:00 2009'

# >>> time.ctime(1256290000)
# 'Fri Oct 23 11:26:40 2009'
                          	  
def compare_results(g1, g2, t0 = 1253700000, delta_t = 2590000, iterations = 20) :
    print "Comparing %d path tree results." % (iterations)
    stations = [v.label for v in g1.vertices if v.label[0:3] == 'sta']
    for i in range(iterations) :
        #origin = stations.pop(random.randrange(len(stations)))
        origin = stations[random.randrange(len(stations))]
        print "Single-source full shortest path tree from %s" % origin
        t = t0 + random.randrange(delta_t)
        print "At %s" % time.ctime(t)
        spt1 = g1.shortest_path_tree( origin,  None, State(1, t))
        spt2 = g2.shortest_path_tree( origin,  None, State(1, t))
        compare_contents(spt1, spt2)
        spt1.destroy()
        spt2.destroy()
        
def test_serialization(gsdb_filename) :
    print "Loading graph database from %s..." % gsdb_filename
    gdb = GraphDatabase(gsdb_filename)
    g1 = gdb.incarnate()
    print "Serializing graph database..."
    g1.serialize('test.gsbin', 'test.mmf')
    print "Duplicating graph from disk..."
    g2 = Graph() # could put gsbin parameter in constructor
    g2.load('test.gsbin', 'test.mmf')
    # do some damage to test that the tests catch errors
    # remove a vertex:
    # g2.remove_vertex(g2.vertices[random.randrange(len(g2.vertices))].label, 1, 1)
    # add an edge:
    # vf = g2.vertices[random.randrange(len(g2.vertices))].label
    # vt = g2.vertices[random.randrange(len(g2.vertices))].label
    # g2.add_edge(vf, vt, Street("walk", 3.5))
    print "Comparing graph contents before/after serialization..."
    compare_contents(g1, g2)
    compare_results(g1, g2)
    print "Cleaning up..."
    g1.destroy()
    g2.destroy()
    
n_discrepancies = 0
n_real_errors = 0
test_serialization('gsdata/bart.linked.gsdb')
test_serialization('gsdata/trimet_13sep2009.linked.gsdb')
test_serialization('gsdata/sfmuni.linked.gsdb')
print "%d real errors in %d discrepancies." % (n_real_errors, n_discrepancies)


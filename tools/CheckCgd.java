/* Native Systre parser audit; no net identification or embedding is requested.
 * javac -cp Systre-19.6.0.jar CheckCgd.java
 * java -cp Systre-19.6.0.jar:. CheckCgd graph.cgd             # unstitched baseline
 * java -cp Systre-19.6.0.jar:. CheckCgd --stitched graph.cgd  # stitched export
 *
 * Baseline mode asserts the 0.3.0 contract: every edge joins consecutive nodes
 * of one segment and every component is finite. Stitched output breaks that by
 * design, so --stitched keeps only the Systre-side checks (parse, declared
 * degrees) and prints Systre's own primitive-cell graph as tab-separated
 * records. tools/check_stitched_cgd.py compares those records to the receipt.
 */
import java.util.*;
import org.gavrog.joss.pgraphs.io.*;
import org.gavrog.joss.pgraphs.basic.*;

public class CheckCgd {
    public static void main(String[] args) throws Exception {
        boolean stitched = args.length == 2 && args[0].equals("--stitched");
        if (!stitched && args.length != 1) throw new IllegalArgumentException("usage: CheckCgd [--stitched] graph.cgd");
        Net net = new NetParser(args[args.length-1]).parseNet();
        if (net == null || !net.isOk()) {
            if (net != null) net.getErrors().forEachRemaining(System.err::println);
            throw new IllegalStateException("CGD parser failed");
        }
        Map<Integer,Integer> degrees = new TreeMap<>();
        for (INode n: net.nodes()) {
            int expected = (Integer) net.getNodeInfo(n,NetParser.CONNECTIVITY);
            if (n.degree() != expected) throw new IllegalStateException("Degree mismatch: "+net.getNodeName(n));
            degrees.merge(n.degree(),1,Integer::sum);
        }
        if (stitched) {
            System.out.println("NET\t"+net.getName()+"\t"+net.getGivenGroup());
            System.out.println("NODES\t"+net.numberOfNodes());
            for (IEdge e: net.edges())
                System.out.println("EDGE\t"+net.getNodeName(e.source())+"\t"+net.getNodeName(e.target())+"\t"+net.getShift(e));
            for (PeriodicGraph.Component c: net.connectedComponents())
                System.out.println("COMPONENT\t"+c.getDimension()+"\t"+c.getGraph().numberOfNodes()+"\t"+c.getMultiplicity());
            net.getWarnings().forEachRemaining(w -> System.out.println("WARNING\t"+w));
            System.out.println("PASS\tparse and declared degrees");
            return;
        }
        for (IEdge e: net.edges()) {
            String a = net.getNodeName(e.source()).split("_N")[0];
            String b = net.getNodeName(e.target()).split("_N")[0];
            if (!a.equals(b)) throw new IllegalStateException("Cross-segment edge: "+e);
            int i = Integer.parseInt(net.getNodeName(e.source()).split("_N")[1]);
            int j = Integer.parseInt(net.getNodeName(e.target()).split("_N")[1]);
            if (Math.abs(i-j)!=1) throw new IllegalStateException("Nonconsecutive trace edge: "+e);
        }
        Map<Integer,Integer> sizes = new TreeMap<>();
        for (PeriodicGraph.Component c: net.connectedComponents()) {
            if (c.getDimension()!=0) throw new IllegalStateException("Unexpected periodic strand");
            sizes.merge(c.getGraph().numberOfNodes(),1,Integer::sum);
        }
        System.out.println("PASS native Systre parser: "+net.getName()+" group="+net.getGivenGroup());
        System.out.println("Primitive-cell nodes="+net.numberOfNodes()+" edges="+net.numberOfEdges());
        System.out.println("Degree histogram="+degrees+" finite component sizes="+sizes);
        net.getWarnings().forEachRemaining(w -> System.out.println("Warning: "+w));
    }
}

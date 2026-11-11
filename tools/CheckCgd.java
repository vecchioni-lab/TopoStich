/* Native Systre parser audit; no net identification or embedding is requested.
 * javac -cp Systre-19.6.0.jar CheckCgd.java
 * java -cp Systre-19.6.0.jar:. CheckCgd graph.cgd
 */
import java.util.*;
import org.gavrog.joss.pgraphs.io.*;
import org.gavrog.joss.pgraphs.basic.*;

public class CheckCgd {
    public static void main(String[] args) throws Exception {
        Net net = new NetParser(args[0]).parseNet();
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

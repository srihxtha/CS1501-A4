import java.io.*;
import java.util.*;

public  class AirlineSystem 
{
  	private String [] cityNames = null; //array of city names 
  	private Digraph G = null; // the digraph
  	private static Scanner scan = null; //scanner 
  	private static final int INFINITY = Integer.MAX_VALUE;
	public static void main(String[] args) throws IOException //main method
	{
	   AirlineSystem airline = new AirlineSystem();
	   scan = new Scanner(System.in);
	   while(true){
	   	// all the switch cases
	    switch(airline.menu()){ // displaying the menu WORKS FULLY
	    	case 1:
	    	  airline.readGraph(); //read graph method WORKS FULLY
	    	  break;
	        case 2:
	          airline.showEverything(); // show the entire list of direct routes, distances, and prices WORKS FULLY
	          break;
	        case 3:
	          airline.kruskalMST();  // show MST - WORKS FULLY 
	          break;
	        case 4:
	          airline.shortestPathMiles(); // show shortest path based on miles // WORKS FULLY
	          break;
	        case 5:
	            airline.shortestPathPrice(); // show shortest path based on price // WORKS FULLY
	           	break;
	        case 6:
	        	airline.shortestPathHops(); // show shortest path based on hops // WORKS FULLY
	        	break;
	        case 7:
	        	airline.flightsUnderCost(); // show trips less than or equal to a price entered by user - WORKS
	        	break;
	        case 8:
	          scan.close();
	          System.exit(0); // quit program 
	          break;
	        default:
	          System.out.println("Incorrect option."); // in case they don't pick a valid number 
	      }
	    }
	  }
	  // printing out the menu for user 
	  private int menu(){
	    System.out.println("*********************************");
	    System.out.println("Welcome to Srihitha Airlines!");
	    System.out.println("1. Read data from a file.");
	    System.out.println("2. Display all direct routes, their distance, and their price.");
	    System.out.println("3. Display a minimum spanning tree for service routes based on distance");
	    System.out.println("4. Compute shortest path based on distance.");
	    System.out.println("5. Compute shortest path based on price.");
	    System.out.println("6. Compute shortest path based on number of hops.");
	    System.out.println("7. Show trips that are less than or equal to a budget that you enter");
	    System.out.println("8. Exit");
	    System.out.println("*********************************");
	    System.out.print("Please choose a menu option (1-8): ");

	    int choice = Integer.parseInt(scan.nextLine());
	    return choice;
	  }
	  // read graph method, user enters the file name and the graph is read in. from lab code 
	  private void readGraph() throws IOException {
	    System.out.println("Please enter graph filename:");
	    String fileName = scan.nextLine();
	    Scanner fileScan = new Scanner(new FileInputStream(fileName));
	    int v = Integer.parseInt(fileScan.nextLine());
	    G = new Digraph(v);

	    cityNames = new String[v];
	    for(int i=0; i<v; i++){
	      cityNames[i] = fileScan.nextLine();
	    }

	    while(fileScan.hasNext()){ //while there is a next line, assign them to these variables. create graph 
	      int from = fileScan.nextInt();
	      int to = fileScan.nextInt();
	      int weight = fileScan.nextInt();
	      double cost = fileScan.nextDouble();
	      G.addEdge(new WeightedDirectedEdge(from-1, to-1, weight, cost));
	      G.addEdge(new WeightedDirectedEdge(to-1, from-1, weight, cost));
	    }
	    fileScan.close();
	    System.out.println("Data imported successfully.");
	    System.out.print("Please press ENTER to continue ...");
	    scan.nextLine();
	  }

	 //show everything method. this method shows the entire list of direct routes, with the distance traveled and the cost
	  //some of this code is taken from the lab 
	 private void showEverything()
	 {
	 
	    if(G == null){ //if the user hasn't put in a graph yet 
      		System.out.println("Please import a graph first (option 1).");
      		System.out.print("Please press ENTER to continue ...");
      		scan.nextLine();
    	} else {
    		System.out.println("-----------------------------------------------------------------------------------");
      		System.out.println("Note that the routes are duplicated, one from each beginning city's point of view");
      		System.out.println("-----------------------------------------------------------------------------------");
      		System.out.println(" ");
      		for (int i = 0; i < G.v; i++) { // printing out the routes, from, to, distance, and cost 
       	 		for (WeightedDirectedEdge e : G.adj(i)) {
          			System.out.print("From " + cityNames[i] + " " + "to " + cityNames[e.to()] + ": " + "(" + e.weight() + " miles)" + " " + "($" + e.cost() + ")" +  " ");
          			System.out.println(" ");
        		}
        	System.out.println();
      		}
      		System.out.print("Please press ENTER to continue ...");
      		scan.nextLine();

    	}
	}

		
  
   //mst method - using kruskal's. A LOT of the code is from the textbook. see below 
  public void kruskalMST()
    {
        
        System.out.println("Minumum Spanning tree");
        System.out.println("--------------------------");
        System.out.println("The edges in the MST based on distance follow: ");

        //calling the Kruskal class 
        KruskalMST showMST = new KruskalMST(G);
        //loop through
        for (WeightedDirectedEdge e : showMST.edges()) {
            System.out.println("From " + cityNames[e.to()] + " to " + cityNames[e.from()] + " - " + e.weight() + " miles");
            System.out.println();
      		
        }
        System.out.print("Please press ENTER to continue ...");
      	scan.nextLine();

    }
	

	
    // shows the shortest path based on the price 
	 private void shortestPathPrice() {
    if(G == null){
      System.out.println("Please import a graph first (option 1).");
      System.out.print("Please press ENTER to continue ...");
      scan.nextLine();
    } else {
      for(int i=0; i<cityNames.length; i++){
        System.out.println(i+1 + ": " + cityNames[i]);
      }
      // basic stuff - ask them what cities they want to travel to and such 
      System.out.print("Please enter source city (1-" + cityNames.length + "): ");
      int source = Integer.parseInt(scan.nextLine());
      System.out.print("Please enter destination city (1-" + cityNames.length + "): ");
      int destination = Integer.parseInt(scan.nextLine());
      source--;
      destination--;
      G.bfs(source); //bfs 
      if(!G.marked[destination]){
        System.out.println("There is no route from " + cityNames[source]
                            + " to " + cityNames[destination]);
      } else {
      	G.dijkstrasCost(source, destination); //Dijkstras 
      	Stack<Integer> stack = new Stack<Integer>(); // using stack 
      	int index = destination;
      	while (index != source) // while youre not going from one city to the same city 
      	{
      		stack.push(index);
      		index = G.edgeTo[index]; //make it an edge 
      		
      	}
      	stack.push(source); //push the beginning city 
       
       //printing 
       System.out.println("Lowest price: " + G.costTo[destination]);
       System.out.println("Path, in order: ");
       int cost = -1;
       int costDest;
       double indivCost = 0;
       LinkedList<WeightedDirectedEdge> source2;
       //pop the stack 
       if (stack.size() > 0)
       {
       	  cost = stack.pop();
       } 
       
       while (stack.size() > 0)
       {
   		//loop through 
       	source2 = G.adj[cost];
       	costDest = stack.pop();
       	for (int i = 0; i < source2.size(); i++)
       	{
       		if (source2.get(i).from() == costDest|| source2.get(i).to() == costDest)
       		{
       		 	indivCost = source2.get(i).cost();
       		}
       	}
       	System.out.print(cityNames[cost] + " " + indivCost + " ");
       	cost = costDest;
     

       }
       if (cost >= 0)
       {
       		System.out.println(cityNames[cost]);
       }

      System.out.print("Please press ENTER to continue ...");
      scan.nextLine();
    }
  }
}
	// show the shortest path in terms of distance
	//code adapted from lab 
	  private void shortestPathMiles() {
    if(G == null){
      System.out.println("Please import a graph first (option 1).");
      System.out.print("Please press ENTER to continue ...");
      scan.nextLine();
    } else {
      for(int i=0; i<cityNames.length; i++){
        System.out.println(i+1 + ": " + cityNames[i]);
      }
      System.out.print("Please enter source city (1-" + cityNames.length + "): ");
      int source = Integer.parseInt(scan.nextLine());
      System.out.print("Please enter destination city (1-" + cityNames.length + "): ");
      int destination = Integer.parseInt(scan.nextLine());
      source--;
      destination--;
      //using bfs 
      G.bfs(source);
      if(!G.marked[destination]){
        System.out.println("There is no route from " + cityNames[source]
                            + " to " + cityNames[destination]);
      } else {
      	//using dijkstras
      	G.dijkstras(source, destination);
      	Stack<Integer> stack = new Stack<Integer>(); //using stack 
      	int index = destination;
      	while (index != source) //while you're not going from one city to the same city 
      	{
      		stack.push(index);
      		index = G.edgeTo[index];
      		
      	}
      	stack.push(source);
       
       // printing 
       System.out.println("Shortest distance: " + G.distTo[destination]);
       System.out.println("Path, in order: ");
       int miles = -1;
       int milesDest;
       int distance = 0;
       LinkedList<WeightedDirectedEdge> source1;
       //using a stack 
       if (stack.size() > 0)
       {
       	  miles = stack.pop();
       } else {
       	 
       }
       
       while (stack.size() > 0)
       {
   
       	source1 = G.adj[miles];
       	milesDest = stack.pop();
       	//loop through 
       	for (int i = 0; i < source1.size(); i++)
       	{
       		if (source1.get(i).from() == milesDest || source1.get(i).to() == milesDest);
       		{
       		 	distance = source1.get(i).weight();
       		}
       	}
       	System.out.print(cityNames[miles] + " " + distance + " ");
       	miles = milesDest;
     

       }
       if (miles >= 0) // as long as there is more to go 
       {
       		System.out.println(cityNames[miles]);
       }
       
      }

      System.out.println(" ");
      System.out.print("Please press ENTER to continue ...");
      scan.nextLine();
    }
  }

  	// shows the shortest. path based on hops. some of this code is adapted from the lab 
	private void shortestPathHops() {
    if(G == null){
      System.out.println("Please import a graph first (option 1).");
      System.out.print("Please press ENTER to continue ...");
      scan.nextLine();
    } else {
      for(int i=0; i<cityNames.length; i++){
        System.out.println(i+1 + ": " + cityNames[i]);
      }
      // sometimes i wonder if you guys read these or just check if there are comments
      // printing basic stuff, asking them what cities 
      System.out.print("Please enter source city (1-" + cityNames.length + "): ");
      int source = Integer.parseInt(scan.nextLine());
      System.out.print("Please enter destination city (1-" + cityNames.length + "): ");
      int destination = Integer.parseInt(scan.nextLine());
      source--;
      destination--;
      G.bfs(source); //bfs 
      if(!G.marked[destination]){
        System.out.println("There is no route from " + cityNames[source]
                            + " to " + cityNames[destination]);
      } else {
      	Stack<Integer> stack = new Stack<Integer>(); //stack 
      	int index = destination;
      	while (index != source) //while youre not going from one city to that same city 
      	{
      		stack.push(index);
      		index = G.edgeTo[index];
      		
      	}
      	stack.push(source);
       //printing stuff 
       System.out.println("Number of hops: " + G.distTo[destination]);
       System.out.println("Path, in order: ");
       while (stack.size() > 0) //while stack size isn't 0
       {
       	System.out.print(cityNames[stack.pop()] + " "); //pop it from stack 
       }
      }
      
      System.out.println(" ");
      System.out.print("Please press ENTER to continue ...");
      scan.nextLine();
    }
  }

  // shows flights that are under a price entered by the user 
  private void flightsUnderCost()
  {
  	LinkedList<WeightedDirectedEdge>[] copyOfadj = G.returnDigraph();
  	Scanner userPrice = new Scanner(System.in);
  	System.out.println("please enter your maximum budget for your trip");

  	double budget = userPrice.nextDouble(); // the user's budget 
  	//printing stuff 
  	System.out.println(" ");
  	System.out.println("Here are all the trips at or under your budget of $" + budget);
  	System.out.println(" ");

  	for (int i = 0; i < G.v; i++)
  	{
  		for (int k = 0; k < copyOfadj[i].size(); k++)
  		{
  			WeightedDirectedEdge route = copyOfadj[i].get(k);
  			int beginCity = route.from(); //city you start in 
  			int destCity = route.to(); // city you end up in 
  			G.dijkstrasCost(beginCity, destCity);
  			if (G.costTo[destCity] <= budget) // as long as the cost isn't above budget
  			{
  				System.out.println("Cost of trip: $" + G.costTo[destCity] + " ");
  				Stack<Integer> stack = new Stack<Integer>();
  				for (int h = destCity; h != beginCity; h = (int)G.edgeTo[h])
  				{
  					stack.push(h);
  				}
  				int displayVertex = beginCity;

  				while (!stack.empty()) // while the stack isn't empty, pop 
  				{
  					int vertPopped = stack.pop();
  					System.out.print(cityNames[beginCity] + " " + (G.costTo[vertPopped] - G.costTo[displayVertex]) + " " + cityNames[vertPopped] + " ");
  					System.out.println(" ");
  					displayVertex = vertPopped;
  				}
  				System.out.println(" ");
  			}
  		}
  	}

  	 System.out.print("Please press ENTER to continue ...");
     scan.nextLine();

  }

// beginnign of classes and helpers. 
  //Digraph comes from lab, added the costTo array though 
  private class Digraph {
    private final int v;
    private int e;
    public LinkedList<WeightedDirectedEdge>[] adj;
    private boolean[] marked;  // marked[v] = is there an s-v path
    private int[] edgeTo;      // edgeTo[v] = previous edge on shortest s-v path
    public int[] distTo;      // distTo[v] = number of edges shortest s-v path
    private double[] costTo;	// costTo array
   // private int[][] graph;


    /**
    * Create an empty digraph with v vertices.
    */
    public Digraph(int v) {
      if (v < 0) throw new RuntimeException("Number of vertices must be nonnegative");
      this.v = v;
      this.e = 0;
      @SuppressWarnings("unchecked")
      LinkedList<WeightedDirectedEdge>[] temp =
      (LinkedList<WeightedDirectedEdge>[]) new LinkedList[v];
      adj = temp;
      for (int i = 0; i < v; i++)
      {
        adj[i] = new LinkedList<WeightedDirectedEdge>();
      }
    }

    /**
    * Add the edge e to this digraph.
    */
    public void addEdge(WeightedDirectedEdge edge) {
      int from = edge.from();
      adj[from].add(edge);
      e++;
    }



    /**
    * Return the edges leaving vertex v as an Iterable.
    * To iterate over the edges leaving vertex v, use foreach notation:
    * <tt>for (WeightedDirectedEdge e : graph.adj(v))</tt>.
    */
    public Iterable<WeightedDirectedEdge> adj(int v) {
      return adj[v];
    }

    public LinkedList<WeightedDirectedEdge>[] returnDigraph()
    {
    	return adj;
    }

    //bfs also comes from lab 
    public void bfs(int source) {
      marked = new boolean[this.v];
      distTo = new int[this.e];
      edgeTo = new int[this.v];

      Queue<Integer> q = new LinkedList<Integer>();
      for (int i = 0; i < v; i++){
        distTo[i] = INFINITY;
        marked[i] = false;
      }
      distTo[source] = 0;
      marked[source] = true;
      q.add(source);

      while (!q.isEmpty()) {
        int v = q.remove();
        for (WeightedDirectedEdge w : adj(v)) {
          if (!marked[w.to()]) {
            edgeTo[w.to()] = v;
            distTo[w.to()] = distTo[v] + 1;
            marked[w.to()] = true;
            q.add(w.to());
          }
        }
      }
    }

    // from lab 
    public void dijkstras(int source, int destination) {
      marked = new boolean[this.v];
      distTo = new int[this.v];
      edgeTo = new int[this.v];


      for (int i = 0; i < v; i++){
        distTo[i] = INFINITY;
        marked[i] = false;
      }
      distTo[source] = 0;
      marked[source] = true;
      int nMarked = 1;

      int current = source;
      while (nMarked < this.v) {
        for (WeightedDirectedEdge w : adj(current)) {
          if (distTo[current]+w.weight() < distTo[w.to()]) {
	      //TODO:update edgeTo and distTo
          	distTo[w.to()] = distTo[current]+w.weight();
          	edgeTo[w.to()] = current;
	      
          }
        }
        //Find the vertex with minimim path distance
        //This can be done more effiently using a priority queue!
        int min = INFINITY;
        current = -1;

        for(int i=0; i<distTo.length; i++){
          if(marked[i])
            continue;
          if(distTo[i] < min){
            min = distTo[i];
            current = i;
          }
        }

	//TODO: Update marked[] and nMarked. Check for disconnected graph.
        if (current == -1)
        {
        	break;
        }
        nMarked++;
        marked[current] = true;
      }
    }

    // made this one - dijkstras but to find the cost of a trip, very similar to the one for distance but with costTo instead 
    public void dijkstrasCost(int source, int destination) {
      marked = new boolean[this.v];
      costTo = new double[this.v];
      edgeTo = new int[this.v];


      for (int i = 0; i < v; i++){
        costTo[i] = INFINITY;
        marked[i] = false;
      }
      costTo[source] = 0;
      marked[source] = true;
      int nMarked = 1;

      int current = source;
      while (nMarked < this.v) {
        for (WeightedDirectedEdge w : adj(current)) {
          if (costTo[current]+w.cost() < costTo[w.to()]) {
	      //TODO:update edgeTo and distTo
          	costTo[w.to()] = costTo[current]+w.cost();
          	edgeTo[w.to()] = current;
	      
          }
        }
        //Find the vertex with minimim path distance
        //This can be done more effiently using a priority queue!
        double min = INFINITY;
        current = -1;

        for(int i=0; i<costTo.length; i++){
          if(marked[i])
            continue;
          if(costTo[i] < min){
            min = costTo[i];
            current = i;
          }
        }

	// Update marked[] and nMarked. Check for disconnected graph.
        if (current == -1)
        {
        	break;
        }
        nMarked++;
        marked[current] = true;
      }
    }

  }

  /**
  *  The <tt>WeightedDirectedEdge</tt> class represents a weighted edge in an directed graph.
  */
  // from lab 
  private class WeightedDirectedEdge {
    private final int v;
    private final int w;
    private int weight;
    private double cost;

    /**
    * Create a directed edge from v to w with given weight.
    */
    // from lab
    public WeightedDirectedEdge(int v, int w, int weight, double cost) {
      this.v = v;
      this.w = w;
      this.weight = weight;
      this.cost = cost;
    }

    public int from(){
      return v;
    }

    public int to(){
      return w;
    }

    public int weight(){
      return weight;
    }
    public double cost()
    {
    	return cost;
    }
  }

  // everything from here is directly from the github from the textbook, but modified a little to adapt to this. 
  // Includes all the dependencies too. 
  // Citation: Sedgewick, Robert, and Kevin Wayne. Algorithms, Fourth Edition. 2016. (github)
  public class KruskalMST {
    private static final double FLOATING_POINT_EPSILON = 1E-12;

    private int weight;                        // weight of MST
    private myQueue<WeightedDirectedEdge> mst = new myQueue<WeightedDirectedEdge>();  // edges in MST

    /**
     * Compute a minimum spanning tree (or forest) of an edge-weighted graph.
     * @param G the edge-weighted graph
     */
    public KruskalMST(Digraph G) {
        // more efficient to build heap by passing array of edges
        MinPQ<WeightedDirectedEdge> pq = new MinPQ<WeightedDirectedEdge>();
        for (int i =0; i < G.v -1; i++)
        {
        	for (WeightedDirectedEdge e : G.adj(i)) 
        	{
            pq.insert(e);
        	}
        }
        

        // run greedy algorithm
        UF uf = new UF(G.v);
        while (!pq.isEmpty() && mst.size() < G.v - 1) {
            WeightedDirectedEdge e = pq.delMin();
            int v = e.from();
            int w = e.to();
            if (uf.find(v) != uf.find(w)) { // v-w does not create a cycle
                uf.union(v, w);  // merge v and w components
                mst.enqueue(e);  // add edge e to mst
                weight += e.weight();
            }
        }

        // check optimality conditions
        assert check(G);
    }

    /**
     * Returns the edges in a minimum spanning tree (or forest).
     * @return the edges in a minimum spanning tree (or forest) as
     *    an iterable of edges
     */
    public Iterable<WeightedDirectedEdge> edges() {
        return mst;
    }

    /**
     * Returns the sum of the edge weights in a minimum spanning tree (or forest).
     * @return the sum of the edge weights in a minimum spanning tree (or forest)
     */
    public double weight() {
        return weight;
    }
    
    // check optimality conditions (takes time proportional to E V lg* V)
    private boolean check(Digraph G) {

        // check total weight
        double total = 0.0;
        for (WeightedDirectedEdge e : edges()) {
            total += e.weight();
        }
        if (Math.abs(total - weight()) > FLOATING_POINT_EPSILON) {
            System.err.printf("Weight of edges does not equal weight(): %f vs. %f\n", total, weight());
            return false;
        }

        // check that it is acyclic
        UF uf = new UF(G.v);
        for (WeightedDirectedEdge e : edges()) {
            int v = e.from(), w = e.to();
            if (uf.find(v) == uf.find(w)) {
                System.err.println("Not a forest");
                return false;
            }
            uf.union(v, w);
        }

        // check that it is a spanning forest
        for (WeightedDirectedEdge e : edges()) {
            int v = e.from(), w = e.to();
            if (uf.find(v) != uf.find(w)) {
                System.err.println("Not a spanning forest");
                return false;
            }
        }

        // check that it is a minimal spanning forest (cut optimality conditions)
        for (WeightedDirectedEdge e : edges()) {

            // all edges in MST except e
            uf = new UF(G.v);
            for (WeightedDirectedEdge f : mst) {
                int x = f.from(), y = f.to();
                if (f != e) uf.union(x, y);
            }
            
            // check that e is min weight edge in crossing cut
            for (WeightedDirectedEdge f : edges()) {
                int x = f.from(), y = f.to();
                if (uf.find(x) != uf.find(y)) {
                    if (f.weight() < e.weight()) {
                        System.err.println("Edge " + f + " violates cut optimality conditions");
                        return false;
                    }
                }
            }

        }

        return true;
    }

  	
}
	public class UF {

    private int[] parent;  // parent[i] = parent of i
    private byte[] rank;   // rank[i] = rank of subtree rooted at i (never more than 31)
    private int count;     // number of components

    /**
     * Initializes an empty union-find data structure with
     * {@code n} elements {@code 0} through {@code n-1}.
     * Initially, each elements is in its own set.
     *
     * @param  n the number of elements
     * @throws IllegalArgumentException if {@code n < 0}
     */
    public UF(int n) {
        if (n < 0) throw new IllegalArgumentException();
        count = n;
        parent = new int[n];
        rank = new byte[n];
        for (int i = 0; i < n; i++) {
            parent[i] = i;
            rank[i] = 0;
        }
    }

    /**
     * Returns the canonical element of the set containing element {@code p}.
     *
     * @param  p an element
     * @return the canonical element of the set containing {@code p}
     * @throws IllegalArgumentException unless {@code 0 <= p < n}
     */
    public int find(int p) {
        validate(p);
        while (p != parent[p]) {
            parent[p] = parent[parent[p]];    // path compression by halving
            p = parent[p];
        }
        return p;
    }

    /**
     * Returns the number of sets.
     *
     * @return the number of sets (between {@code 1} and {@code n})
     */
    public int count() {
        return count;
    }
  
    /**
     * Returns true if the two elements are in the same set.
     *
     * @param  p one element
     * @param  q the other element
     * @return {@code true} if {@code p} and {@code q} are in the same set;
     *         {@code false} otherwise
     * @throws IllegalArgumentException unless
     *         both {@code 0 <= p < n} and {@code 0 <= q < n}
     * @deprecated Replace with two calls to {@link #find(int)}.
     */
    @Deprecated
    public boolean connected(int p, int q) {
        return find(p) == find(q);
    }
  
    /**
     * Merges the set containing element {@code p} with the 
     * the set containing element {@code q}.
     *
     * @param  p one element
     * @param  q the other element
     * @throws IllegalArgumentException unless
     *         both {@code 0 <= p < n} and {@code 0 <= q < n}
     */
    public void union(int p, int q) {
        int rootP = find(p);
        int rootQ = find(q);
        if (rootP == rootQ) return;

        // make root of smaller rank point to root of larger rank
        if      (rank[rootP] < rank[rootQ]) parent[rootP] = rootQ;
        else if (rank[rootP] > rank[rootQ]) parent[rootQ] = rootP;
        else {
            parent[rootQ] = rootP;
            rank[rootP]++;
        }
        count--;
    }

    // validate that p is a valid index
    private void validate(int p) {
        int n = parent.length;
        if (p < 0 || p >= n) {
            throw new IllegalArgumentException("index " + p + " is not between 0 and " + (n-1));  
        }
    }
    }

    public class MinPQ<Key> implements Iterable<Key> {
    private Key[] pq;                    // store items at indices 1 to n
    private int n;                       // number of items on priority queue
    private Comparator<Key> comparator;  // optional comparator

    /**
     * Initializes an empty priority queue with the given initial capacity.
     *
     * @param  initCapacity the initial capacity of this priority queue
     */
    public MinPQ(int initCapacity) {
        pq = (Key[]) new Object[initCapacity + 1];
        n = 0;
    }

    /**
     * Initializes an empty priority queue.
     */
    public MinPQ() {
        this(1);
    }

    /**
     * Initializes an empty priority queue with the given initial capacity,
     * using the given comparator.
     *
     * @param  initCapacity the initial capacity of this priority queue
     * @param  comparator the order in which to compare the keys
     */
    public MinPQ(int initCapacity, Comparator<Key> comparator) {
        this.comparator = comparator;
        pq = (Key[]) new Object[initCapacity + 1];
        n = 0;
    }

    /**
     * Initializes an empty priority queue using the given comparator.
     *
     * @param  comparator the order in which to compare the keys
     */
    public MinPQ(Comparator<Key> comparator) {
        this(1, comparator);
    }

    /**
     * Initializes a priority queue from the array of keys.
     * <p>
     * Takes time proportional to the number of keys, using sink-based heap construction.
     *
     * @param  keys the array of keys
     */
    public MinPQ(Key[] keys) {
        n = keys.length;
        pq = (Key[]) new Object[keys.length + 1];
        for (int i = 0; i < n; i++)
            pq[i+1] = keys[i];
        for (int k = n/2; k >= 1; k--)
            sink(k);
        assert isMinHeap();
    }

    /**
     * Returns true if this priority queue is empty.
     *
     * @return {@code true} if this priority queue is empty;
     *         {@code false} otherwise
     */
    public boolean isEmpty() {
        return n == 0;
    }

    /**
     * Returns the number of keys on this priority queue.
     *
     * @return the number of keys on this priority queue
     */
    public int size() {
        return n;
    }

    /**
     * Returns a smallest key on this priority queue.
     *
     * @return a smallest key on this priority queue
     * @throws NoSuchElementException if this priority queue is empty
     */
    public Key min() {
        if (isEmpty()) throw new NoSuchElementException("Priority queue underflow");
        return pq[1];
    }

    // resize the underlying array to have the given capacity
    private void resize(int capacity) {
        assert capacity > n;
        Key[] temp = (Key[]) new Object[capacity];
        for (int i = 1; i <= n; i++) {
            temp[i] = pq[i];
        }
        pq = temp;
    }

    /**
     * Adds a new key to this priority queue.
     *
     * @param  x the key to add to this priority queue
     */
    public void insert(Key x) {
        // double size of array if necessary
        if (n == pq.length - 1) resize(2 * pq.length);

        // add x, and percolate it up to maintain heap invariant
        pq[++n] = x;
        swim(n);
        assert isMinHeap();
    }

    /**
     * Removes and returns a smallest key on this priority queue.
     *
     * @return a smallest key on this priority queue
     * @throws NoSuchElementException if this priority queue is empty
     */
    public Key delMin() {
        if (isEmpty()) throw new NoSuchElementException("Priority queue underflow");
        Key min = pq[1];
        exch(1, n--);
        sink(1);
        pq[n+1] = null;     // to avoid loitering and help with garbage collection
        if ((n > 0) && (n == (pq.length - 1) / 4)) resize(pq.length / 2);
        assert isMinHeap();
        return min;
    }


   /***************************************************************************
    * Helper functions to restore the heap invariant.
    ***************************************************************************/

    private void swim(int k) {
        while (k > 1 && greater(k/2, k)) {
            exch(k, k/2);
            k = k/2;
        }
    }

    private void sink(int k) {
        while (2*k <= n) {
            int j = 2*k;
            if (j < n && greater(j, j+1)) j++;
            if (!greater(k, j)) break;
            exch(k, j);
            k = j;
        }
    }

   /***************************************************************************
    * Helper functions for compares and swaps.
    ***************************************************************************/
    private boolean greater(int i, int j) {
        
    	if (((WeightedDirectedEdge) pq[i]).weight() > ((WeightedDirectedEdge)pq[j]).weight())
    	{
    		return true;
    	} else {
    		return false;
    	}
        
    }

    private void exch(int i, int j) {
        Key swap = pq[i];
        pq[i] = pq[j];
        pq[j] = swap;
    }

    // is pq[1..n] a min heap?
    private boolean isMinHeap() {
        for (int i = 1; i <= n; i++) {
            if (pq[i] == null) return false;
        }
        for (int i = n+1; i < pq.length; i++) {
            if (pq[i] != null) return false;
        }
        if (pq[0] != null) return false;
        return isMinHeapOrdered(1);
    }

    // is subtree of pq[1..n] rooted at k a min heap?
    private boolean isMinHeapOrdered(int k) {
        if (k > n) return true;
        int left = 2*k;
        int right = 2*k + 1;
        if (left  <= n && greater(k, left))  return false;
        if (right <= n && greater(k, right)) return false;
        return isMinHeapOrdered(left) && isMinHeapOrdered(right);
    }


    /**
     * Returns an iterator that iterates over the keys on this priority queue
     * in ascending order.
     * <p>
     * The iterator doesn't implement {@code remove()} since it's optional.
     *
     * @return an iterator that iterates over the keys in ascending order
     */
    public Iterator<Key> iterator() {
        return new HeapIterator();
    }

    private class HeapIterator implements Iterator<Key> {
        // create a new pq
        private MinPQ<Key> copy;

        // add all items to copy of heap
        // takes linear time since already in heap order so no keys move
        public HeapIterator() {
            if (comparator == null) copy = new MinPQ<Key>(size());
            else                    copy = new MinPQ<Key>(size(), comparator);
            for (int i = 1; i <= n; i++)
                copy.insert(pq[i]);
        }

        public boolean hasNext()  { return !copy.isEmpty();                     }
        public void remove()      { throw new UnsupportedOperationException();  }

        public Key next() {
            if (!hasNext()) throw new NoSuchElementException();
            return copy.delMin();
        }
    }
}
    public class myQueue<Item> implements Iterable<Item> {
    private Node<Item> first;    // beginning of queue
    private Node<Item> last;     // end of queue
    private int n;               // number of elements on queue

    // helper linked list class
    private class Node<Item> {
        private Item item;
        private Node<Item> next;
    }

    /**
     * Initializes an empty queue.
     */
    public myQueue() {
        first = null;
        last  = null;
        n = 0;
    }

    /**
     * Returns true if this queue is empty.
     *
     * @return {@code true} if this queue is empty; {@code false} otherwise
     */
    public boolean isEmpty() {
        return first == null;
    }

    /**
     * Returns the number of items in this queue.
     *
     * @return the number of items in this queue
     */
    public int size() {
        return n;
    }

    /**
     * Returns the item least recently added to this queue.
     *
     * @return the item least recently added to this queue
     * @throws NoSuchElementException if this queue is empty
     */
    public Item peek() {
        if (isEmpty()) throw new NoSuchElementException("Queue underflow");
        return first.item;
    }

    /**
     * Adds the item to this queue.
     *
     * @param  item the item to add
     */
    public void enqueue(Item item) {
        Node<Item> oldlast = last;
        last = new Node<Item>();
        last.item = item;
        last.next = null;
        if (isEmpty()) first = last;
        else           oldlast.next = last;
        n++;
    }

    /**
     * Removes and returns the item on this queue that was least recently added.
     *
     * @return the item on this queue that was least recently added
     * @throws NoSuchElementException if this queue is empty
     */
    public Item dequeue() {
        if (isEmpty()) throw new NoSuchElementException("Queue underflow");
        Item item = first.item;
        first = first.next;
        n--;
        if (isEmpty()) last = null;   // to avoid loitering
        return item;
    }

    /**
     * Returns a string representation of this queue.
     *
     * @return the sequence of items in FIFO order, separated by spaces
     */
    public String toString() {
        StringBuilder s = new StringBuilder();
        for (Item item : this) {
            s.append(item);
            s.append(' ');
        }
        return s.toString();
    } 

    /**
     * Returns an iterator that iterates over the items in this queue in FIFO order.
     *
     * @return an iterator that iterates over the items in this queue in FIFO order
     */
    public Iterator<Item> iterator()  {
        return new LinkedIterator(first);  
    }

    // an iterator, doesn't implement remove() since it's optional
    private class LinkedIterator implements Iterator<Item> {
        private Node<Item> current;

        public LinkedIterator(Node<Item> first) {
            current = first;
        }

        public boolean hasNext()  { return current != null;                     }
        public void remove()      { throw new UnsupportedOperationException();  }

        public Item next() {
            if (!hasNext()) throw new NoSuchElementException();
            Item item = current.item;
            current = current.next; 
            return item;
        }
    }

}
}
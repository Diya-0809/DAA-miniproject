# DAA Analysis — Algorithmic Relationship Mapping of a Family Tree

**Project:** Algorithmic Relationship Mapping via Graph Data Structures
**Implementation:** `family_tree.cpp` (C++14, standard library only)
**Model:** Directed Acyclic Graph (DAG) stored as a bidirectional N-ary adjacency list

---

## 1. Project Overview & Data Structure Selection

### 1.1 Project Concept

The project represents a family lineage as a directed graph and answers
genealogical questions by composing four classical graph algorithms:

| Query | Algorithm | Implementation entry point |
|---|---|---|
| Kinship classification (siblings, cousins, aunt/uncle, …) | Lowest Common Ancestor (LCA) on a DAG | `findLca()` → `relationshipText()` |
| All people exactly $k$ generations below a person (e.g. grandchildren) | Breadth-First Search (BFS) downwards | `generationsBelow()` |
| Every unbroken bloodline back to an origin ancestor | Depth-First Search (DFS) upwards | `lineagePaths()` |
| Degree of separation / consanguinity between two relatives | Unweighted shortest path (BFS on the undirected view) | `shortestKinshipPath()` → `degreesOfSeparation()` |

Throughout this document the standard graph-theoretic notation is used:

- $V$ = the **vertices** of the graph, i.e. the *people* (records containing an
  identifier, a name, an optional photo path, and adjacency lists);
- $E$ = the **edges** of the graph, i.e. the *parent–child relationships*, each
  directed from parent to child.

Every edge points exactly one generation downwards, so the graph can contain no
directed cycle; moreover each person has **in-degree $\leq 2$** (at most two
parents), which implies

$$|E| \leq 2|V| \quad\Longrightarrow\quad O(V + E) = O(V),$$

i.e. the graph is sparse and every $O(V+E)$ traversal is effectively linear in
the number of people. Two auxiliary parameters appear in the analysis of the
depth-first traversal: $d$ (the pedigree depth, i.e. the number of generations
under consideration) and $P$ (the number of distinct root-to-query bloodlines,
which is an *output-size* parameter and may be exponential — see §2.3).

### 1.2 Why a Directed Acyclic Graph?

A family tree is conventionally drawn as a tree, but it is **not** a tree in the
formal sense. Three structural facts force the DAG model:

1. **Two parents per child.** Every person has (at most) two parents, so every
   vertex has in-degree two. A tree — and a binary tree in particular — requires
   a single unique parent for every non-root vertex; the moment a child has a
   mother *and* a father, the parent-pointer structure branches upwards and the
   "one parent" invariant is violated.
2. **Merging lineages (multiple paths).** Because of (1), two people may reach a
   shared ancestor through *several distinct paths*. For example, full siblings
   share both parents, and their grandchildren are reachable from a
   great-grandparent through four different chains. A structure in which a vertex
   can be reached by more than one path from a common ancestor is by definition a
   DAG, not a tree. This is precisely why the *lowest* common ancestor is
   non-trivial and must be computed algorithmically rather than by walking two
   unique root paths.
3. **Acyclicity.** Every edge points strictly from a higher generation to the
   next lower one. Since generations are well ordered, no directed cycle can
   exist; the graph is therefore acyclic *by construction*, and the edge-insertion
   routine additionally verifies this invariant with a reachability test
   (`validateEdge()`, $O(V+E)$) so that malformed input can never corrupt it.

Directedness matters as well: the direction *parent → child* encodes the
biological meaning of the relationship. The same edge traversed downwards means
"descendant", traversed upwards means "ancestor"; reversing it would silently
swap the two roles. Keeping the direction explicit lets one algorithm exploit it
(LCA, lineage DFS) while another deliberately symmetrises it (degree of
separation, §2.4).

### 1.3 Why an Adjacency List Rather Than an Array or a Matrix?

The graph is stored as an **N-ary adjacency list in two directions**:

```
Person {
    id, name, photo,
    parents[]   // reverse adjacency: child  -> parents   (in-degree  ≤ 2)
    children[]  // forward adjacency: parent -> children  (out-degree unbounded)
}
```

The alternative representations were rejected for the following reasons.

**Compared with a flat array (or a hash-indexed collection of records):**

- An array supports only *positional* access. Finding "the grandchildren of person
  $x$" or "the closest shared ancestor of $x$ and $y$" requires scanning the whole
  array and re-deriving the relationships from scratch, i.e. $\Theta(V^2)$ work if
  every record is compared against every other record, or repeated $\Theta(V)$
  scans per intermediate step.
- An array cannot *share* structure: recording both parents of a child means
  duplicating ancestry blocks, and any update (adding a spouse, correcting a
  link) forces data movement. In the adjacency list, adding an edge is $O(1)$
  amortised (two `push_back`s).
- An array has no notion of *distance*, so concepts such as "generation" or
  "degree of separation" have no representation at all, whereas on an adjacency
  list they fall out of BFS levels for free.

**Compared with a binary tree:**

- As argued in §1.2(1), a binary tree models at most one parent per child. Two
  parents per child require an N-ary (general) graph.
- Even if the graph *were* a tree, a binary tree imposes an artificial left/right
  partition on siblings: a person with five children would need dummy internal
  nodes purely to satisfy the degree-2 out-degree constraint, inflating $V$
  without adding information.
- Tree algorithms assume **unique paths**. In a family DAG a vertex may be
  reachable from an ancestor by several paths of *different* lengths (a person
  may be simultaneously a nephew and a grand-nephew of the same ancestor), so
  naïve tree formulas for LCA or depth return the wrong answer. The algorithms
  must instead reason about *shortest* generational distance, which is exactly
  what BFS provides.
- Query sets differ: a binary search tree answers membership queries in
  $O(h)$, but a family tree is not ordered by key — it is traversed by
  *relationship*, which is a graph reachability question.

**Compared with an adjacency matrix:**

- The graph is sparse ($|E| \leq 2|V|$, average degree $\leq 2$). An adjacency
  matrix requires $\Theta(V^2)$ space — for a genealogical database of
  $10^6$ people that is $10^{12}$ entries, versus $\Theta(V+E) = \Theta(V)$
  entries for the adjacency list.
- Scanning the neighbours of a vertex in a matrix costs $\Theta(V)$, so a full
  traversal costs $\Theta(V^2)$; in the adjacency list it costs
  $\Theta(V+E)$.

**Summary of the selection:**

| Criterion | Array / record list | Binary tree | Adjacency matrix | **DAG + N-ary adjacency list (chosen)** |
|---|---|---|---|---|
| Space | $O(V)$ | $O(V)$ | $O(V^2)$ | $O(V+E) = O(V)$ |
| Enumerate neighbours of a person | $O(V)$ | $O(1)$–$O(V)$ | $O(V)$ | $O(\deg)$ |
| Full traversal (BFS/DFS) | $O(V^2)$ | $O(V)$ | $O(V^2)$ | $O(V+E)$ |
| Two parents per child | duplicated data | impossible | possible | **native (in-degree ≤ 2)** |
| Multiple paths to an ancestor | n/a | forbidden | possible | **modelled, queried via LCA** |
| Generational distance | not represented | depth only | computable | **BFS levels, $O(V+E)$** |

---

## 2. Algorithmic Analysis (The Core Report)

Notation used in this section: $V$ = vertices (people), $E$ = edges
(parent–child ties), $h$ = height of the pedigree, $d$ = depth of the explored
region, $P$ = number of distinct bloodlines ending at the queried person.
Recall $|E| \le 2|V|$, so every bound stated as $O(V+E)$ is in practice
$O(V)$.

### 2.1 Lowest Common Ancestor (LCA)

#### 2.1.1 Application to the family tree

Given two people $A$ and $B$, their *common ancestors* are the vertices that
occur on some ancestral path of $A$ **and** on some ancestral path of $B$. The
**lowest** such ancestor is the one closest to the pair; it is the vertex through
which their kinship is "bottlenecked", and therefore the natural justification
for the relationship label.

Because the structure is a DAG, several common ancestors may exist at different
depths, so the classical tree algorithm (climb both root paths until the paths
diverge) does not apply. The implementation generalises it as follows:

1. Run an **upward BFS** from $A$ over the reverse adjacency list, obtaining
   $d_A[x]$ = number of generations from $A$ up to each ancestor $x$
   (first discovery = shortest distance on an unweighted graph);
2. run the same upward BFS from $B$, obtaining $d_B[x]$;
3. **intersect** the two maps: every key present in both is a common ancestor;
4. select the common ancestor minimising the total generational span
   $d_A[x] + d_B[x]$, breaking ties by the imbalance $|d_A[x]-d_B[x]|$ (the most
   symmetric choice yields the intuitive answer) and finally by vertex id
   (determinism).

The winning pair $(d_A, d_B)$ is a *relationship fingerprint* that is decoded
into a kinship term:

| $(d_A, d_B)$ | Interpretation | Example output |
|---|---|---|
| $(0, 0)$ | same person | — |
| $(0, k)$, $k \ge 1$ | $A$ is an ancestor of $B$ | "grandparent of …" |
| $(k, 0)$, $k \ge 1$ | $A$ is a descendant of $B$ | "great-grandchild of …" |
| $(1, 1)$ | siblings (in-degree intersection distinguishes **full** vs **half**) | "full siblings" |
| $(1, k)$, $k \ge 2$ | $A$ is collaterally above $B$ | "aunt or uncle of …" |
| $(k, 1)$, $k \ge 2$ | mirror image | "niece or nephew of …" |
| $(k, l)$, $k,l \ge 2$ | collateral relatives: degree $\min(k,l)-1$, removal $|k-l|$ | "$1^{\text{st}}$ cousins, once removed" |

Selecting the ancestor that minimises $d_A + d_B$ is exactly the definition of
"lowest" in a DAG where a unique path does not exist: it is the common ancestor
with the smallest combined generational distance, and the tie-break on
$|d_A - d_B|$ prevents degenerate choices such as labelling two first cousins
through a great-great-grandparent.

#### 2.1.2 Complexity

- **Time — $O(V + E)$.** The dominant work is two upward BFS traversals, each
  $O(V + E)$ in the worst case (every vertex and edge is enqueued at most once).
  The intersection scan is $O(|\text{anc}(A)|) \le O(V)$, and the fingerprint
  decode is $O(1)$ (the parent-set intersection compares at most $2 \times 2$
  entries). Overall: $T_{\text{LCA}}(V,E) = 2 \cdot O(V+E) + O(V) + O(1) = O(V+E)$.
  In practice only the *ancestor cones* of $A$ and $B$ are explored, so the
  realised cost is $O(\text{anc}(A) + \text{anc}(B) + E_A + E_B)$, which is
  substantially smaller than $O(V+E)$ for distant relatives in a large database.
- **Space — $O(V)$.** Two distance maps plus two BFS queues, each bounded by
  $O(V)$; the result record is $O(1)$ apart from the optional list of equally
  good common ancestors ($O(V)$ in the worst case).

*Worst case for a true tree* — if the input is a tree, the two ancestor chains
have length $O(h)$ and the algorithm degrades gracefully to $O(h)$ per chain;
the $O(V+E)$ bound is the honest worst case for the general DAG, where ancestor
cones may cover the whole graph.

### 2.2 Breadth-First Search (BFS)

#### 2.2.1 Application to the family tree

BFS answers **generational tier queries**: *"who are George's grandchildren?"*,
*"list everyone exactly $k$ generations below person $x$."*

The procedure runs a standard forward BFS (parent → child) from the query vertex
`root`. Since every edge carries unit weight (one generation), the BFS level at
which a vertex is **first discovered** equals its shortest generational distance
below `root`. Three details make the result genealogically correct:

- **Shortest distance wins.** In a DAG a person may be reachable at several
  depths (simultaneously a nephew at depth 2 and a grand-nephew at depth 3); the
  first discovery reports the closest one, so each person is counted exactly once
  at their *nearest* tier.
- **Early termination.** Vertices already at level $k$ are not expanded (their
  children would lie at level $k+1$), so the search stops as soon as tier $k$ is
  exhausted rather than traversing the whole graph.
- **Deterministic output.** Vertices at level $k$ are collected and sorted by id.

BFS is the correct choice — rather than DFS — because the query is about
*distance from a single source*, and BFS is the canonical algorithm for shortest
paths in unweighted graphs; a DFS could reach tier $k$ vertices only after
wandering through unrelated branches and would have to explore the entire
subtree to be sure it had found them all.

#### 2.2.2 Complexity

- **Time — $O(V + E)$ worst case.** Each vertex is enqueued at most once and
  each edge is inspected at most once, giving $\sum_{v} O(\deg^+(v)) = O(E)$ for
  the expansion plus $O(V)$ for queue and map operations. With early termination
  at level $k$ the realised cost is $O(V_k + E_k)$, where $V_k, E_k$ are the
  vertices and edges of the explored region (the *ball* of radius $k$ around
  `root`), and the final collection of tier $k$ costs $O(V_k + \log V_k)$
  because of the sort (the scan itself is $O(V_k)$).
- **Space — $O(V)$.** The level map and the queue together hold at most $O(V)$
  entries; $O(V_k)$ when the search is truncated at level $k$.

### 2.3 Depth-First Search (DFS)

#### 2.3.1 Application to the family tree

DFS answers **lineage queries**: *"trace Penelope's bloodlines back to the oldest
known ancestors."* A bloodline is a maximal ancestral path
$$\text{Penelope} \leftarrow \text{father} \leftarrow \text{grandfather} \leftarrow \cdots \leftarrow
\text{origin ancestor},$$
i.e. a path that terminates at a vertex with no recorded parents.

The implementation is a recursive DFS over the **reverse adjacency list**
(child → parents) with an explicit `current` path vector:

- `push_back(u)` on entry,
- if `parents(u)` is empty, the current path is copied into the result set — a
  complete bloodline has been reconstructed,
- otherwise recurse into **every** parent,
- `pop_back()` on exit (**backtracking**), which restores the invariant that
  `current` always contains exactly the vertices on the current search path.

DFS is the natural algorithm here for three reasons: (i) the object of interest
*is* a path, and DFS explores one path to its end before branching, which matches
the recursive structure of ancestry; (ii) DFS needs only the current path plus
the recursion stack, $O(d)$ space, whereas BFS would have to hold an entire
generational frontier, $O(\text{width})$; (iii) because the graph is acyclic,
the recursion always terminates — there is no need for grey/black colouring to
detect back-edges, although a visited set would be required if cycles were
possible.

#### 2.3.2 Complexity

- **Time — output-sensitive: $O\!\left(\sum_{i=1}^{P} |\pi_i|\right) = O(P \cdot d)$.**
  The DFS performs $O(1)$ work per vertex *per path* on which that vertex
  appears; $P$ is the number of distinct bloodlines ending at the queried
  person and $d$ their length. Two regimes must be distinguished:
  - *Reporting one arbitrary lineage* (or "is there a line to ancestor $x$"):
    a single DFS without path enumeration, $O(V + E)$.
  - *Enumerating all bloodlines* (the implemented query): the bound is
    $O(P \cdot d)$ with $P$ potentially **exponential** in the pedigree depth.
    If every ancestor has two parents, then $P$ can reach $2^{d}$, so
    $T = O(d \cdot 2^{d})$. This is not an artefact of the implementation:
    counting the root-to-node paths in a DAG is a #P-complete problem, so the
    output itself is exponentially large and any algorithm must spend at least
    time linear in its output ($\Omega(P \cdot d)$ just to print the paths).
    Memoising completed subtrees reduces repeated *computation* of shared
    ancestry to $O(V+E)$, but the materialised list of paths remains
    $\Theta(P \cdot d)$ — the blow-up lives in the data, not in the code.
- **Space — $O(d)$ for the traversal, $O(P \cdot d)$ for the result.**
  The recursion depth and the `current` path vector are each bounded by
  $d \le V$ (the graph is acyclic, so a path never repeats a vertex), and the
  returned set of bloodlines stores $P$ paths of length up to $d$. Intermediate
  space is therefore $O(V)$ in the worst case — markedly better than the
  $O(V)$-sized frontier BFS would need for a wide generation.

### 2.4 Shortest Path (Unweighted / BFS)

#### 2.4.1 Application to the family tree

This query computes the **degree of separation** (informally, consanguinity):
*"how many links are there in the shortest chain of family ties between these
two distant cousins?"*

Kinship ties are logically bidirectional — "my mother's brother" and "my uncle"
are the same statement — so for this query every directed edge is treated as
undirected, giving the symmetrised neighbourhood
$$N(u) = \text{parents}(u) \cup \text{children}(u).$$
A plain BFS is then run from $A$ in this undirected view, storing for every
discovered vertex its distance and its predecessor. When $B$ is dequeued, the
predecessor chain is walked back to $A$ and reversed, returning the actual chain
(e.g. *Penelope → Leo → George → Arthur*) together with its length
$\ell = |\text{path}| - 1$, the degree of separation. If $B$ is never reached,
the two vertices lie in different connected components and the query reports
$-1$.

**Why BFS and not Dijkstra or Bellman–Ford?** All edges have unit weight, so
BFS's greedy level-by-level expansion *is* Dijkstra with a FIFO queue: vertices
are visited in non-decreasing distance order, and the first time $B$ is seen its
distance is minimal. The specialised algorithms would only add overhead —
$O(E + V\log V)$ with a binary heap (or $O(VE)$ for Bellman–Ford) versus
$O(V+E)$ for BFS — without improving correctness.

**Correctness sketch.** By induction on the BFS rounds: after round $r$ every
vertex at distance $r$ (and no vertex at distance $> r$) has been discovered,
because any shortest $A \leadsto v$ path has a prefix of length $r$ ending at a
vertex discovered in round $r-1$. Hence the predecessor chain reconstructed for
$B$ is a shortest path, and its edge count is the minimum number of family links
joining $A$ and $B$.

#### 2.4.2 Complexity

- **Time — $O(V + E)$ worst case.** Each vertex is enqueued once and each of its
  at most four neighbours ($\le 2$ parents, $\le 2$ in the general out-degree
  case) is examined once, so the inner loop is
  $\sum_{v} O(\deg^+(v) + \deg^-(v)) = O(E)$, plus $O(V)$ for map operations and
  $O(\ell)$ to rebuild the path, where $\ell \le V-1$ is the path length.
  Because the search stops as soon as $B$ is dequeued, the realised cost is
  $O(V_{\text{ball}} + E_{\text{ball}})$ over the explored neighbourhood —
  $O(1)$ when the two people are close relatives, $O(V+E)$ when they are in
  different, weakly connected branches.
- **Space — $O(V)$.** Distance map, predecessor map and queue each hold $O(V)$
  entries; the returned path adds $O(\ell) \subseteq O(V)$.

#### 2.4.3 Comparative complexity of the four queries

| Algorithm | Application | Time | Space |
|---|---|---|---|
| LCA (two upward BFS + intersection) | shared ancestors → relationship type | $O(V + E)$ | $O(V)$ |
| BFS downwards | generational tier (grandchildren, $k$-th descendants) | $O(V + E)$ worst case; $O(V_k + E_k)$ with early stop | $O(V)$ |
| DFS upwards | enumerate all bloodlines to origin ancestors | $O(P \cdot d)$ output-sensitive; $O(V+E)$ for one lineage | $O(d)$ stack + $O(P \cdot d)$ output |
| Shortest path (BFS, undirected view) | degree of separation / consanguinity | $O(V + E)$ | $O(V)$ |

*(Edge insertion is guarded by an $O(V+E)$ reachability check so the acyclicity
of the DAG is maintained at all times; the guard itself never exceeds the cost
of the queries it protects.)*

#### 2.4.4 Empirical validation

The implementation includes a self-test (`--test`, 19 assertions covering all
five routines — all passing) and a benchmark (`--bench`) that builds random
family DAGs in which person $i$ receives 1–2 parents from $\{0,\dots,i-1\}$
(a topological order therefore exists by construction, so acyclicity is
guaranteed) and measures the average cost of 2000 random queries per size:

| $V$ (people) | build (ms) | LCA ($\mu$s) | degrees of separation ($\mu$s) | BFS downwards $k{=}3$ ($\mu$s) |
|-----------:|-----------:|---------:|-------------:|-------------------:|
| 500 | 0.00 | 11.77 | 186.85 | 1.17 |
| 1000 | 0.00 | 31.52 | 343.80 | 5.72 |
| 2000 | 21.37 | 11.05 | 683.69 | 11.41 |
| 4000 | 24.75 | 61.74 | 2776.92 | 7.73 |
| 8000 | 89.97 | 92.01 | 7214.46 | 13.84 |

The measurements are consistent with the asymptotic analysis. The
degree-of-separation query grows roughly linearly with $V$ (random pairs are
usually far apart, so BFS explores most of the connected component), the
tier query with early termination stays one to two orders of magnitude cheaper
than a full traversal, and the LCA cost varies with the *size of the ancestor
cones* of the particular random pair rather than monotonically with $V$ — exactly
the instance-sensitive behaviour predicted by the $O(\text{anc}(A)+\text{anc}(B))$
refinement in §2.1.2.

---

## 3. Conclusion

The four algorithms examined above form a coherent toolkit in which each
traversal is matched to the shape of the question it answers: **BFS** exploits
its level structure to measure *generational distance* — both downwards (tier
queries such as "all grandchildren") and upwards/undirected (shortest chains for
degrees of separation) — guaranteeing that the first visit to any vertex is its
closest; **DFS** exploits its depth-first, backtracking character to *reconstruct
complete paths*, enumerating every bloodline back to an origin ancestor in
$O(d)$ working space and reporting honestly that the enumeration is
output-sensitive ($O(P\cdot d)$, exponential only because the answer itself can
be); and the **LCA** composition of two upward BFS passes converts raw reachability
into a *semantic* relationship label (sibling, aunt/uncle, $n^{\text{th}}$
cousin $k$ times removed) in $O(V+E)$ time and $O(V)$ space.

Two design decisions make these bounds attainable in practice. First, the
**DAG + bidirectional adjacency list** representation: with $|E| \le 2|V|$ the
graph is sparse, so every $O(V+E)$ traversal is effectively linear in the number
of people, and storing both `parents[]` and `children[]` lets upward and
downward queries share the same $O(\deg)$ neighbour access instead of paying a
search. Second, the **invariant-preserving edge guard**, which keeps the graph
acyclic at $O(V+E)$, ensuring the correctness arguments above (shortest-distance
= first discovery, guaranteed termination of DFS) remain valid for the lifetime
of the database.

For real-world genealogical databases — which routinely reach $10^{6}$–$10^{9}$
records — these properties translate directly into interactive query latency:
kinship, tier and separation queries are single-pass linear scans over a
sub-graph, memory footprint is $\Theta(V+E)$ rather than $\Theta(V^2)$, and no
query ever requires the quadratic all-pairs comparisons that a naïve record-list
implementation would need. Where the pure analysis becomes the bottleneck — the
exponential enumeration of *all* bloodlines (§2.3) — the remedy is algorithmic
rather than structural: memoised DP over the DAG, early termination for
existential variants of the query, or precomputed ancestor indices, all of which
preserve the $O(V+E)$ per-query framework established here. The result is a
model in which the classical graph algorithms of the DAA curriculum are not mere
exercises but the load-bearing components of a scalable relationship-mapping
system.

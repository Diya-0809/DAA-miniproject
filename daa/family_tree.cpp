/**
 * ============================================================================
 *  Algorithmic Relationship Mapping — A Family Tree as a Directed Acyclic
 *  Graph (DAG)
 * ============================================================================
 *  Course : Design & Analysis of Algorithms (DAA)
 *  Language : C++14, standard library only, single translation unit
 *
 *  Build : g++ -std=c++14 -O2 -Wall -Wextra family_tree.cpp -o family_tree
 *  Run   : ./family_tree           interactive menu
 *          ./family_tree --test    validation demo + sample queries (self-check)
 *          ./family_tree --bench   timing experiment (V vs. microseconds)
 *          ./family_tree --help    usage
 *
 *  GRAPH MODEL
 *  -----------
 *    Vertex V : a Person (id, name, photo path, parents, children)
 *    Edge   E : a directed parent -> child tie. Every edge points exactly one
 *               generation downwards, so the graph can never contain a
 *               directed cycle. Because a person may have up to two parents,
 *               the structure is in general a DAG, not a tree: a descendant
 *               can be reached through two different lines of ancestors.
 *
 *  STORAGE — N-ARY ADJACENCY LIST (two directions)
 *  ------------------------------------------------
 *    children[] : forward adjacency (parent -> child) — BFS downwards,
 *                 reachability test for cycle rejection.
 *    parents[]  : reverse adjacency (child -> parent) — upward BFS for the
 *                 LCA, DFS lineage tracing, and the neighbourhood union used
 *                 by the shortest-path query.
 *
 *  ALGORITHMS IMPLEMENTED (see DAA Analysis.md for the full analysis)
 *  -----------------------------------------------------------------
 *    [1] Lowest Common Ancestor on a DAG  -> relationship classification
 *    [2] Breadth-First Search downwards   -> generational tier queries
 *    [3] Depth-First Search upwards       -> bloodline / lineage paths
 *    [4] Unweighted shortest path (BFS)   -> degrees of separation
 *    [5] Edge validation (keeps the graph a legal family DAG)
 * ============================================================================
 */

#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// ============================================================================
//  Person — the vertex payload
// ============================================================================
struct Person {
    int id = -1;
    std::string name;
    std::string photo;            // e.g. "photos/arthur.jpg" (may be empty)
    std::vector<int> parents;     // in-degree  <= 2 (enforced by validation)
    std::vector<int> children;    // out-degree unbounded (N-ary)
};

// ============================================================================
//  FamilyTree — owns the DAG and every algorithm
// ============================================================================
class FamilyTree {
public:
    // ------------------------------------------------------------------------
    //  Vertex insertion. TIME O(1) amortised, SPACE O(1) for the new record.
    // ------------------------------------------------------------------------
    bool addPerson(int id, const std::string& name, const std::string& photo = "") {
        if (members_.count(id)) return false;
        Person p;
        p.id = id;
        p.name = name;
        p.photo = photo;
        members_[id] = p;
        return true;
    }

    bool hasPerson(int id) const { return members_.count(id) != 0; }
    std::size_t size() const { return members_.size(); }
    const Person& person(int id) const { return members_.at(id); }

    // ========================================================================
    //  ALGORITHM 5 — EDGE VALIDATION (preserves the family-DAG invariant)
    // -----------------------------------------------------------------------
    //  Before inserting parent -> child we guarantee:
    //    (a) both endpoints exist,          (b) no self-loop,
    //    (c) no duplicate edge,             (d) in-degree of the child < 2,
    //    (e) no cycle: adding parent->child creates a cycle iff `parent` is
    //        already reachable from `child` by walking DOWN the graph.
    //
    //  TIME  O(V + E) worst case — check (e) is a reachability BFS over the
    //        forward adjacency list; (a)-(d) are O(1) because in-degree <= 2.
    //  SPACE O(V) — visited set and queue of that BFS.
    // ========================================================================
    enum class Verdict { Ok, NoSuchPerson, SelfLoop, DuplicateEdge,
                         ParentLimit, Cycle };

    static std::string verdictText(Verdict v) {
        switch (v) {
            case Verdict::Ok:           return "edge accepted";
            case Verdict::NoSuchPerson: return "REJECTED: no such person";
            case Verdict::SelfLoop:     return "REJECTED: a person cannot be their own parent";
            case Verdict::DuplicateEdge:return "REJECTED: edge already exists";
            case Verdict::ParentLimit:  return "REJECTED: a person may have at most 2 parents";
            case Verdict::Cycle:        return "REJECTED: edge would create a cycle";
            default:                    return "REJECTED: unknown error";
        }
    }

    Verdict validateEdge(int parentId, int childId) const {
        if (!members_.count(parentId) || !members_.count(childId))
            return Verdict::NoSuchPerson;
        if (parentId == childId)
            return Verdict::SelfLoop;

        const Person& child = members_.at(childId);
        if (std::find(child.parents.begin(), child.parents.end(), parentId)
                != child.parents.end())
            return Verdict::DuplicateEdge;
        if (child.parents.size() >= kMaxParents)
            return Verdict::ParentLimit;

        // (e) reachability test: BFS from the child DOWN to the parent.
        std::unordered_set<int> visited;
        std::queue<int> frontier;
        visited.insert(childId);
        frontier.push(childId);
        while (!frontier.empty()) {
            int u = frontier.front();
            frontier.pop();
            if (u == parentId) return Verdict::Cycle;
            for (int v : members_.at(u).children)
                if (visited.insert(v).second) frontier.push(v);
        }
        return Verdict::Ok;
    }

    bool addParentChild(int parentId, int childId, Verdict& verdict) {
        verdict = validateEdge(parentId, childId);
        if (verdict != Verdict::Ok) return false;
        members_[childId].parents.push_back(parentId);
        members_[parentId].children.push_back(childId);
        return true;
    }

    // ========================================================================
    //  UPWARD BFS — building block of the LCA algorithm
    // -----------------------------------------------------------------------
    //  BFS from `start` against the direction of the edges (child -> parent).
    //  Each upward edge costs exactly one generation, so the BFS level at
    //  which an ancestor is FIRST discovered equals its generational distance
    //  from `start` — i.e. BFS yields shortest distances on an unweighted
    //  graph. The vertex `start` itself receives distance 0.
    //
    //  TIME  O(V + E) worst case; in practice only the ancestor cone of
    //        `start` is explored: O(A + E_A) where A is the number of its
    //        ancestors and E_A the number of edges among them.
    //  SPACE O(V) — distance map + queue.
    // ========================================================================
    std::unordered_map<int, int> ancestorDepths(int start) const {
        std::unordered_map<int, int> depth;
        if (!members_.count(start)) return depth;
        std::queue<int> q;
        depth[start] = 0;
        q.push(start);
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (int p : members_.at(u).parents) {
                if (depth.insert(std::make_pair(p, depth[u] + 1)).second)
                    q.push(p);
            }
        }
        return depth;
    }

    // ========================================================================
    //  ALGORITHM 1 — LOWEST COMMON ANCESTOR (LCA) ON A DAG
    // -----------------------------------------------------------------------
    //  In a tree the LCA is unique and computable by walking up both root
    //  paths; in a DAG (two parents per child) a pair may share SEVERAL
    //  common ancestors at different depths, so we:
    //    1. run upward BFS from A  ->  da[x] = generations from A up to x,
    //    2. run upward BFS from B  ->  db[x] = generations from B up to x,
    //    3. intersect the two maps  -> every common ancestor of A and B,
    //    4. keep the common ancestor minimising da + db (closest in total
    //       generations); ties are broken by |da - db| (most symmetric, which
    //       gives the intuitive answer such as "grandparent" for two cousins)
    //       and finally by the smallest id, for determinism.
    //  The winning pair (da, db) is the relationship fingerprint that
    //  relationshipText() decodes (siblings, cousins N-times removed, ...).
    //
    //  TIME  O(V + E): two upward BFS traversals plus one linear scan of the
    //        smaller distance map. The traversals dominate.
    //  SPACE O(V): the two distance maps and the two BFS queues.
    // ========================================================================
    struct LcaInfo {
        int ancestor = -1;     // chosen LCA id; -1 = no common ancestor
        int da = 0;            // generations from A up to the LCA
        int db = 0;            // generations from B up to the LCA
        int shared = 0;        // number of common ancestors at the optimal depth
        bool found() const { return ancestor != -1; }
    };

    LcaInfo findLca(int a, int b) const {
        std::unordered_map<int, int> da = ancestorDepths(a);
        std::unordered_map<int, int> db = ancestorDepths(b);

        LcaInfo best;
        int bestScore = INT_MAX, bestGap = INT_MAX;

        for (std::unordered_map<int, int>::const_iterator it = da.begin();
             it != da.end(); ++it) {
            std::unordered_map<int, int>::const_iterator jt = db.find(it->first);
            if (jt == db.end()) continue;               // not an ancestor of B

            int score = it->second + jt->second;        // total generations up
            int gap   = std::abs(it->second - jt->second);
            bool better =
                score < bestScore ||
                (score == bestScore && gap < bestGap) ||
                (score == bestScore && gap == bestGap && it->first < best.ancestor);
            if (better) {
                best.ancestor = it->first;
                best.da = it->second;
                best.db = jt->second;
                bestScore = score;
                bestGap = gap;
            }
        }

        // Count every common ancestor lying at the optimal depth; this is what
        // separates FULL siblings (two shared parents) from HALF siblings.
        if (best.found()) {
            for (std::unordered_map<int, int>::const_iterator it = da.begin();
                 it != da.end(); ++it) {
                std::unordered_map<int, int>::const_iterator jt = db.find(it->first);
                if (jt != db.end() && it->second + jt->second == bestScore)
                    ++best.shared;
            }
        }
        return best;
    }

    // ------------------------------------------------------------------------
    //  Linguistic helpers: decode the (da, db) fingerprint into a kinship term
    //  TIME O(1) — pure string construction; parent-set intersection is O(4).
    // ------------------------------------------------------------------------
    static std::string ordinal(int n) {
        std::string s = std::to_string(n);
        int h = n % 100;
        if (h >= 11 && h <= 13) return s + "th";
        switch (n % 10) {
            case 1:  return s + "st";
            case 2:  return s + "nd";
            case 3:  return s + "rd";
            default: return s + "th";
        }
    }
    static std::string greats(int g) {
        std::string s;
        for (int i = 0; i < g; ++i) s += "great-";
        return s;
    }
    static std::string ancestorTerm(int d) {          // d generations up
        if (d == 1) return "parent";
        if (d == 2) return "grandparent";
        return greats(d - 2) + "grandparent";
    }
    static std::string descendantTerm(int d) {        // d generations down
        if (d == 1) return "child";
        if (d == 2) return "grandchild";
        return greats(d - 2) + "grandchild";
    }
    static std::string avuncularTerm(int d) {         // aunt / uncle (d >= 2)
        if (d == 2) return "aunt or uncle";
        return greats(d - 2) + "aunt or uncle";
    }
    static std::string collateralTerm(int d) {        // niece / nephew
        if (d == 2) return "niece or nephew";
        return greats(d - 2) + "niece or nephew";
    }

    std::string describe(int a, int b, const LcaInfo& lca) const {
        const std::string& na = members_.at(a).name;
        const std::string& nb = members_.at(b).name;
        int da = lca.da, db = lca.db;

        if (da == 0 && db == 0)
            return na + " and " + nb + " are the same person";

        if (da == 0)   // A is an ancestor of B
            return na + " is the " + ancestorTerm(db) + " of " + nb;

        if (db == 0)   // A is a descendant of B
            return na + " is the " + descendantTerm(da) + " of " + nb;

        if (da == 1 && db == 1) {
            // Siblings — how many parents do they share?
            const std::vector<int>& pa = members_.at(a).parents;
            const std::vector<int>& pb = members_.at(b).parents;
            int common = 0;
            for (std::size_t i = 0; i < pa.size(); ++i)
                if (std::find(pb.begin(), pb.end(), pa[i]) != pb.end()) ++common;
            if (common >= 2) return na + " and " + nb + " are full siblings";
            return na + " and " + nb + " are half siblings";
        }

        if (da == 1 && db >= 2)          // A is one tier below the LCA
            return na + " is the " + avuncularTerm(db) + " of " + nb;

        if (da >= 2 && db == 1)          // B is one tier below the LCA
            return na + " is the " + collateralTerm(da) + " of " + nb;

        // Both >= 2 : collateral relatives (cousins).
        //   degree  = min(da, db) - 1,  removed = |da - db|
        int degree = std::min(da, db) - 1;
        int removed = std::abs(da - db);
        std::string core = na + " and " + nb + " are " + ordinal(degree) +
                           " cousins";
        if (removed == 0) return core;
        std::string times;
        if (removed == 1)      times = "once";
        else if (removed == 2) times = "twice";
        else                   times = std::to_string(removed) + " times";
        return core + ", " + times + " removed";
    }

    // ------------------------------------------------------------------------
    //  Public relationship query: LCA + naming + the ancestor that justified it
    //  TIME O(V + E)   SPACE O(V)   (inherited from findLca)
    // ------------------------------------------------------------------------
    std::string relationshipText(int a, int b) const {
        if (!members_.count(a) || !members_.count(b))
            return "ERROR: unknown person id.";
        LcaInfo lca = findLca(a, b);
        if (!lca.found())
            return members_.at(a).name + " and " + members_.at(b).name +
                   " share no common ancestor (different families)";

        std::ostringstream out;
        out << describe(a, b, lca);
        out << "\n    closest common ancestor: " << members_.at(lca.ancestor).name
            << " (" << lca.ancestor << ")  [d(" << a << ")=" << lca.da
            << ", d(" << b << ")=" << lca.db << ", shared at this depth = "
            << lca.shared << "]";
        return out.str();
    }

    // ========================================================================
    //  ALGORITHM 2 — BREADTH-FIRST SEARCH DOWNWARDS (generational tier query)
    // -----------------------------------------------------------------------
    //  Forward BFS (parent -> child) from `root`. Because all edges are
    //  unweighted, the level at which a vertex is FIRST discovered equals its
    //  shortest generational distance below `root`. In a DAG a person may be
    //  reachable at several depths (e.g. both nephew and grand-nephew); the
    //  BFS reports the closest one only, so each vertex is counted once.
    //  The search stops expanding at depth k, then collects level k.
    //
    //  Example: generationsBelow(George, 2) -> all grandchildren of George.
    //
    //  TIME  O(V + E) worst case (k >= height); O(V_k + E_k) when the tier is
    //        reached early, since expansion stops at depth k.
    //  SPACE O(V) — level map + queue; O(V_k) in the early-stopped case.
    // ========================================================================
    std::vector<int> generationsBelow(int root, int k) const {
        std::vector<int> tier;
        if (k <= 0 || !members_.count(root)) return tier;

        std::unordered_map<int, int> level;
        std::queue<int> q;
        level[root] = 0;
        q.push(root);

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            if (level[u] >= k) continue;              // children would exceed k
            for (std::size_t i = 0; i < members_.at(u).children.size(); ++i) {
                int c = members_.at(u).children[i];
                if (level.insert(std::make_pair(c, level[u] + 1)).second)
                    q.push(c);
            }
        }
        for (std::unordered_map<int, int>::const_iterator it = level.begin();
             it != level.end(); ++it)
            if (it->second == k) tier.push_back(it->first);
        std::sort(tier.begin(), tier.end());
        return tier;
    }

    // ========================================================================
    //  ALGORITHM 3 — DEPTH-FIRST SEARCH UPWARDS (bloodline / lineage tracing)
    // -----------------------------------------------------------------------
    //  Recursive DFS with backtracking over the REVERSE adjacency list
    //  (child -> parents). Every root-to-start path that terminates at a
    //  person with no recorded parents is one complete, unbroken bloodline
    //  back to an origin ancestor; it is copied into the result set. The
    //  current path is maintained with push_back / pop_back.
    //
    //  TIME  Let P = number of distinct bloodlines (root paths) ending at
    //        `start`. The DFS visits each vertex on each such path once, so
    //        T = O(sum over paths of path length) = O(P * d) where d is the
    //        pedigree depth; with memoisation of completed subtrees this
    //        drops to O(V + E). P is the output size and CAN BE
    //        EXPONENTIAL: if every ancestor has 2 parents and the pedigree
    //        has depth d, P can reach 2^d — the classic exponential blow-up
    //        of path enumeration on a DAG (counting paths is #P-complete).
    //        Reporting a single arbitrary lineage instead of all of them is
    //        O(V + E) — one DFS without backtracking enumeration.
    //  SPACE O(d) recursion depth + O(d) for the current path (d <= V),
    //        plus O(P * d) for the materialised result.
    // ========================================================================
    std::vector<std::vector<int> > lineagePaths(int start) const {
        std::vector<std::vector<int> > paths;
        if (!members_.count(start)) return paths;
        std::vector<int> current;

        std::function<void(int)> dfs = [&](int u) {
            current.push_back(u);
            const std::vector<int>& pars = members_.at(u).parents;
            if (pars.empty()) {
                paths.push_back(current);             // origin ancestor reached
            } else {
                for (std::size_t i = 0; i < pars.size(); ++i) dfs(pars[i]);
            }
            current.pop_back();                       // backtrack
        };
        dfs(start);
        return paths;
    }

    // Pretty-print: "Penelope <- Leo <- George <- Arthur"
    std::string lineageToString(const std::vector<int>& path) const {
        std::ostringstream out;
        for (std::size_t i = 0; i < path.size(); ++i) {
            if (i) out << "  <-  ";
            out << members_.at(path[i]).name << " (" << path[i] << ")";
        }
        return out.str();
    }

    // ========================================================================
    //  ALGORITHM 4 — UNWEIGHTED SHORTEST PATH (BFS on the undirected view)
    // -----------------------------------------------------------------------
    //  Degree of separation / consanguinity ignores the direction of a union:
    //  kinship ties work both ways ("my mother's brother" = "my uncle").
    //  We therefore run a plain BFS on the symmetrised neighbourhood
    //      N(u) = parents(u) U children(u)
    //  and stop when the target is dequeued. BFS on an unweighted graph
    //  visits vertices in non-decreasing distance order, so the first time the
    //  target is seen its distance is minimal. A predecessor map rebuilds the
    //  actual chain of ties as the returned path.
    //
    //  Returns an empty path (degree -1) if the two people are in different
    //  connected components.
    //
    //  TIME  O(V + E) worst case; the search stops as soon as the target is
    //        reached, so it is O(V_b + E_b) over the explored ball.
    //  SPACE O(V) — distance/predecessor map + queue (+ O(l) for a path of
    //        length l when reconstructing).
    // ========================================================================
    bool shortestKinshipPath(int a, int b, std::vector<int>& path) const {
        path.clear();
        if (!members_.count(a) || !members_.count(b)) return false;
        if (a == b) { path.push_back(a); return true; }

        std::unordered_map<int, int> dist;
        std::unordered_map<int, int> prev;
        std::queue<int> q;
        dist[a] = 0;
        q.push(a);

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            if (u == b) break;
            int du = dist[u];
            int neighbourhood[4];                     // <= 2 parents + <= 2 used
            int n = 0;
            for (std::size_t i = 0; i < members_.at(u).parents.size(); ++i)
                neighbourhood[n++] = members_.at(u).parents[i];
            for (std::size_t i = 0; i < members_.at(u).children.size(); ++i)
                neighbourhood[n++] = members_.at(u).children[i];
            for (int i = 0; i < n; ++i) {
                int v = neighbourhood[i];
                if (dist.insert(std::make_pair(v, du + 1)).second) {
                    prev[v] = u;
                    q.push(v);
                }
            }
        }

        if (!dist.count(b)) return false;             // disconnected
        for (int v = b; v != a; v = prev[v]) path.push_back(v);
        path.push_back(a);
        std::reverse(path.begin(), path.end());
        return true;
    }

    // Degree of separation = number of edges in the shortest chain = |path| - 1
    int degreesOfSeparation(int a, int b) const {
        std::vector<int> path;
        if (!shortestKinshipPath(a, b, path)) return -1;
        return static_cast<int>(path.size()) - 1;
    }

    // ========================================================================
    //  GRAPHVIZ .dot EXPORT  (visualisation of the DAG)
    //  TIME O(V + E)   SPACE O(V)   (one forward BFS assigns generations)
    // ========================================================================
    bool exportDot(const std::string& file) const {
        std::unordered_map<int, int> depth;
        std::queue<int> q;
        for (std::unordered_map<int, Person>::const_iterator it = members_.begin();
             it != members_.end(); ++it)
            if (it->second.parents.empty()) {
                depth[it->first] = 0;
                q.push(it->first);
            }
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (std::size_t i = 0; i < members_.at(u).children.size(); ++i) {
                int c = members_.at(u).children[i];
                if (depth.insert(std::make_pair(c, depth[u] + 1)).second) q.push(c);
            }
        }

        static const char* palette[] = {"#ffd9b3", "#b3e6ff", "#c9f7c9",
                                        "#e6ccff", "#ffcccc"};
        std::ofstream out(file.c_str());
        if (!out) return false;

        out << "// Family tree as a DAG (parent -> child)\n"
            << "// Render: dot -Tpng " << file << " -o family.png\n"
            << "digraph FamilyTree {\n"
            << "  rankdir=TB;\n"
            << "  node [shape=box, style=\"rounded,filled\", "
               "fontname=\"Helvetica\"];\n"
            << "  edge [arrowhead=none];\n\n";

        std::vector<int> ids = sortedIds();
        for (std::size_t i = 0; i < ids.size(); ++i) {
            const Person& p = members_.at(ids[i]);
            int d = depth.count(ids[i]) ? depth[ids[i]] : 0;
            out << "  n" << ids[i] << " [label=\"" << escapeDot(p.name)
                << "\\n(" << ids[i] << ")\", fillcolor=\""
                << palette[d % 5] << "\"";
            if (!p.photo.empty())
                out << ", image=\"" << escapeDot(p.photo)
                    << "\", labelloc=\"b\"";
            out << "];\n";
        }
        out << "\n";
        for (std::size_t i = 0; i < ids.size(); ++i)
            for (std::size_t j = 0; j < members_.at(ids[i]).children.size(); ++j)
                out << "  n" << ids[i] << " -> n"
                    << members_.at(ids[i]).children[j] << ";\n";
        out << "}\n";
        return true;
    }

    std::vector<int> sortedIds() const {
        std::vector<int> ids;
        for (std::unordered_map<int, Person>::const_iterator it = members_.begin();
             it != members_.end(); ++it)
            ids.push_back(it->first);
        std::sort(ids.begin(), ids.end());
        return ids;
    }

    static std::string escapeDot(const std::string& s) {
        std::string r;
        for (std::size_t i = 0; i < s.size(); ++i) {
            if (s[i] == '"' || s[i] == '\\') r += '\\';
            r += s[i];
        }
        return r;
    }

private:
    static const int kMaxParents = 2;
    std::unordered_map<int, Person> members_;   // adjacency-list vertex store
};

// ============================================================================
//  SAMPLE 4-GENERATION FAMILY
//  ------------------------------------------------------------------
//  >>> Replace with your own family: addPerson(...) + addParentChild(...) <<<
//
//  Generation 1 (ids 1-6)   : three couples (the great-grandparents)
//  Generation 2 (ids 7-11)  : their children (grandparents)
//  Generation 3 (ids 12-16) : cousins' parents; id 14 is a HALF-sibling
//  Generation 4 (ids 17-19) : the youngest generation
// ============================================================================
FamilyTree buildSampleFamily() {
    FamilyTree t;

    // Generation 1
    t.addPerson(1, "Arthur Anderson",   "photos/arthur_anderson.jpg");
    t.addPerson(2, "Beatrice Anderson", "photos/beatrice_anderson.jpg");
    t.addPerson(3, "Charles Miller",    "photos/charles_miller.jpg");
    t.addPerson(4, "Dorothy Miller",    "photos/dorothy_miller.jpg");
    t.addPerson(5, "Edward Smith",      "photos/edward_smith.jpg");
    t.addPerson(6, "Frances Smith",     "photos/frances_smith.jpg");
    // Generation 2
    t.addPerson(7,  "George Anderson",  "photos/george_anderson.jpg");
    t.addPerson(8,  "Helen Anderson",   "photos/helen_anderson.jpg");
    t.addPerson(9,  "Irene Miller",     "photos/irene_miller.jpg");
    t.addPerson(10, "Jack Smith",       "photos/jack_smith.jpg");
    t.addPerson(11, "Karen Smith",      "photos/karen_smith.jpg");
    // Generation 3
    t.addPerson(12, "Leo Anderson",     "photos/leo_anderson.jpg");
    t.addPerson(13, "Mia Anderson",     "photos/mia_anderson.jpg");
    t.addPerson(14, "Oscar Anderson",   "photos/oscar_anderson.jpg"); // half-sib
    t.addPerson(15, "Quentin Reed",     "photos/quentin_reed.jpg");
    t.addPerson(16, "Rachel Reed",      "photos/rachel_reed.jpg");
    // Generation 4
    t.addPerson(17, "Penelope Anderson","photos/penelope_anderson.jpg");
    t.addPerson(18, "Samuel Anderson",  "photos/samuel_anderson.jpg");
    t.addPerson(19, "Tara Anderson",    "photos/tara_anderson.jpg");

    FamilyTree::Verdict v;
    int edges[][2] = {
        {1, 7},  {2, 7},   {1, 8},  {2, 8},           // Arthur+Beatrice -> George, Helen
        {3, 9},  {4, 9},                              // Charles+Dorothy -> Irene
        {5, 10}, {6, 10}, {5, 11}, {6, 11},           // Edward+Frances  -> Jack, Karen
        {7, 12}, {9, 12}, {7, 13}, {9, 13},           // George+Irene    -> Leo, Mia
        {7, 14},                                      // George          -> Oscar (1 parent)
        {10, 15}, {11, 15}, {10, 16}, {11, 16},       // Jack+Karen      -> Quentin, Rachel
        {12, 17}, {15, 17},                           // Leo+Quentin     -> Penelope
        {13, 18}, {16, 18},                           // Mia+Rachel      -> Samuel
        {14, 19}, {16, 19}                            // Oscar+Rachel    -> Tara
    };
    int nEdges = static_cast<int>(sizeof(edges) / sizeof(edges[0]));
    for (int i = 0; i < nEdges; ++i)
        if (!t.addParentChild(edges[i][0], edges[i][1], v))
            std::cerr << "sample edge rejected: " << edges[i][0] << "->"
                      << edges[i][1] << " : " << FamilyTree::verdictText(v)
                      << "\n";
    return t;
}

// ============================================================================
//  SELF-TEST  (./family_tree --test)
//  Runs the validation demo and the sample queries, checks the expected
//  answers, and exits with status 1 if any check fails.
// ============================================================================
namespace {
int gChecks = 0;
int gFailures = 0;

void check(const std::string& label, const std::string& got,
           const std::string& want) {
    ++gChecks;
    bool ok = (got == want);
    if (!ok) ++gFailures;
    std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << label << "\n"
              << "         got : " << got << "\n";
    if (!ok) std::cout << "         want: " << want << "\n";
}

void checkInt(const std::string& label, long long got, long long want) {
    check(label, std::to_string(got), std::to_string(want));
}
}  // namespace

void runSelfTest(FamilyTree& t) {
    std::cout << "==========================================================\n"
              << "  PART A — EDGE VALIDATION (algorithm 5)\n"
              << "==========================================================\n";
    FamilyTree::Verdict v;

    t.addParentChild(14, 12, v);   // third parent for Leo
    check("third parent for Leo (14->12) rejected",
          FamilyTree::verdictText(v), "REJECTED: a person may have at most 2 parents");
    t.addParentChild(12, 1, v);    // cycle: Leo parent of Arthur
    check("cycle 12->1 rejected", FamilyTree::verdictText(v),
          "REJECTED: edge would create a cycle");
    t.addParentChild(7, 7, v);     // self loop
    check("self loop 7->7 rejected", FamilyTree::verdictText(v),
          "REJECTED: a person cannot be their own parent");
    t.addParentChild(7, 12, v);    // duplicate
    check("duplicate 7->12 rejected", FamilyTree::verdictText(v),
          "REJECTED: edge already exists");
    t.addParentChild(99, 12, v);   // unknown person
    check("unknown person 99->12 rejected", FamilyTree::verdictText(v),
          "REJECTED: no such person");

    std::cout << "\n==========================================================\n"
              << "  PART B — ALGORITHM 1: LCA / RELATIONSHIP QUERIES\n"
              << "==========================================================\n";

    // Helpers that extract only the kinship phrase (first line) of the report.
    auto phrase = [&t](int a, int b) {
        std::string s = t.relationshipText(a, b);
        std::size_t nl = s.find('\n');
        return nl == std::string::npos ? s : s.substr(0, nl);
    };

    check("relationship(1, 12) Arthur-Leo",
          phrase(1, 12), "Arthur Anderson is the grandparent of Leo Anderson");
    check("relationship(12, 13) Leo-Mia",
          phrase(12, 13), "Leo Anderson and Mia Anderson are full siblings");
    check("relationship(14, 12) Oscar-Leo",
          phrase(14, 12), "Oscar Anderson and Leo Anderson are half siblings");
    check("relationship(17, 18) Penelope-Samuel",
          phrase(17, 18), "Penelope Anderson and Samuel Anderson are 1st cousins");
    check("relationship(19, 1) Tara-Arthur",
          phrase(19, 1), "Tara Anderson is the great-grandchild of Arthur Anderson");
    check("relationship(17, 9) Penelope-Irene",
          phrase(17, 9), "Penelope Anderson is the grandchild of Irene Miller");

    std::cout << "\n==========================================================\n"
              << "  PART C — ALGORITHM 2: BREADTH-FIRST SEARCH (tiers)\n"
              << "==========================================================\n";

    std::vector<int> tier = t.generationsBelow(1, 2);
    std::ostringstream oss;
    for (std::size_t i = 0; i < tier.size(); ++i) {
        if (i) oss << ", ";
        oss << tier[i];
    }
    check("generationsBelow(Arthur=1, k=2)", oss.str(), "12, 13, 14");

    tier = t.generationsBelow(7, 2);
    oss.str("");
    oss.clear();
    for (std::size_t i = 0; i < tier.size(); ++i) {
        if (i) oss << ", ";
        oss << tier[i];
    }
    check("generationsBelow(George=7, k=2) = his grandchildren",
          oss.str(), "17, 18, 19");

    std::cout << "\n==========================================================\n"
              << "  PART D — ALGORITHM 3: DEPTH-FIRST SEARCH (bloodlines)\n"
              << "==========================================================\n";

    std::vector<std::vector<int> > lines = t.lineagePaths(17);
    checkInt("number of bloodlines of Penelope (17)",
             static_cast<long long>(lines.size()), 8);
    for (std::size_t i = 0; i < lines.size(); ++i)
        std::cout << "    line " << (i + 1) << ": "
                  << t.lineageToString(lines[i]) << "\n";
    bool allEndAtRoot = true;
    for (std::size_t i = 0; i < lines.size(); ++i)
        if (!t.person(lines[i].back()).parents.empty()) allEndAtRoot = false;
    checkInt("every bloodline ends at a parentless origin ancestor",
             allEndAtRoot ? 1 : 0, 1);

    std::cout << "\n==========================================================\n"
              << "  PART E — ALGORITHM 4: SHORTEST PATH (degrees of sep.)\n"
              << "==========================================================\n";

    checkInt("degreesOfSeparation(Penelope=17, Samuel=18)",
             t.degreesOfSeparation(17, 18), 4);
    checkInt("degreesOfSeparation(Penelope=17, Arthur=1)",
             t.degreesOfSeparation(17, 1), 3);
    checkInt("degreesOfSeparation(Penelope=17, Penelope=17)",
             t.degreesOfSeparation(17, 17), 0);

    std::vector<int> path;
    if (t.shortestKinshipPath(17, 1, path)) {
        std::cout << "    shortest chain 17 -> 1: ";
        for (std::size_t i = 0; i < path.size(); ++i) {
            if (i) std::cout << " - ";
            std::cout << t.person(path[i]).name;
        }
        std::cout << "\n";
    }

    std::cout << "\n==========================================================\n"
              << "  PART F — VISUALISATION\n"
              << "==========================================================\n";
    bool written = t.exportDot("family.dot");
    checkInt("exportDot(\"family.dot\")", written ? 1 : 0, 1);
    std::cout << "    render with: dot -Tpng family.dot -o family.png\n";

    std::cout << "\n==========================================================\n"
              << "  " << (gChecks - gFailures) << "/" << gChecks
              << " checks passed";
    if (gFailures) std::cout << "  (" << gFailures << " FAILED)";
    std::cout << "\n==========================================================\n";
}

// ============================================================================
//  BENCHMARK  (./family_tree --bench)
//  Random family DAGs of increasing size; person i receives 1-2 parents from
//  {0..i-1}, so a topological order exists by construction (ids) and the
//  graph is acyclic by construction. Reports average microseconds/query.
// ============================================================================
FamilyTree randomFamily(int n, unsigned seed) {
    FamilyTree t;
    std::mt19937 rng(seed);
    for (int i = 0; i < n; ++i)
        t.addPerson(i, "P" + std::to_string(i));
    FamilyTree::Verdict v;
    for (int i = 1; i < n; ++i) {
        int inDegree = 1 + static_cast<int>(rng() % 2);   // 1 or 2 parents
        for (int k = 0; k < inDegree; ++k) {
            int p = static_cast<int>(rng() % i);          // parent id < child id
            t.addParentChild(p, i, v);                    // duplicates ignored
        }
    }
    return t;
}

void runBench() {
    const int SIZES[] = {500, 1000, 2000, 4000, 8000};
    const int QUERIES = 2000;
    const unsigned SEED = 12345;

    std::cout << "Timing experiment: random family DAGs, " << QUERIES
              << " random queries per cell (average microseconds/query).\n"
              << "Acyclicity guaranteed: every parent id < child id.\n\n";
    std::cout << std::left
              << "| V (people) | build (ms) | LCA (us) | degrees (us) |"
                 " BFS-down k=3 (us) |\n"
              << "|-----------:|-----------:|---------:|-------------:"
                 "|-------------------:|\n";

    for (std::size_t si = 0; si < sizeof(SIZES) / sizeof(SIZES[0]); ++si) {
        int n = SIZES[si];
        std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
        FamilyTree t = randomFamily(n, SEED);
        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        double buildMs =
            std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
                .count() / 1000.0;

        std::mt19937 rng(SEED + n);
        long long lcaUs = 0, degUs = 0, downUs = 0;
        for (int q = 0; q < QUERIES; ++q) {
            int a = static_cast<int>(rng() % n);
            int b = static_cast<int>(rng() % n);
            int c = static_cast<int>(rng() % n);

            std::chrono::steady_clock::time_point s0 =
                std::chrono::steady_clock::now();
            volatile bool sink1 = t.findLca(a, b).found();
            std::chrono::steady_clock::time_point s1 =
                std::chrono::steady_clock::now();
            lcaUs += std::chrono::duration_cast<std::chrono::nanoseconds>(s1 - s0)
                         .count();

            s0 = std::chrono::steady_clock::now();
            volatile int sink2 = t.degreesOfSeparation(a, b);
            s1 = std::chrono::steady_clock::now();
            degUs += std::chrono::duration_cast<std::chrono::nanoseconds>(s1 - s0)
                         .count();

            s0 = std::chrono::steady_clock::now();
            std::vector<int> sink3 = t.generationsBelow(c, 3);
            s1 = std::chrono::steady_clock::now();
            downUs += std::chrono::duration_cast<std::chrono::nanoseconds>(s1 - s0)
                          .count();
            (void)sink1;
            (void)sink2;
            (void)sink3;
        }

        std::cout << "| " << std::right << std::setw(11) << n
                  << " | " << std::setw(10) << std::fixed
                  << std::setprecision(2) << buildMs
                  << " | " << std::setw(8) << (lcaUs / 1000.0 / QUERIES)
                  << " | " << std::setw(13) << (degUs / 1000.0 / QUERIES)
                  << " | " << std::setw(19)
                  << (downUs / 1000.0 / QUERIES) << " |\n";
    }
}

// ============================================================================
//  INTERACTIVE MENU
// ============================================================================
void printMenu() {
    std::cout <<
        "\n-------- ALGORITHMIC RELATIONSHIP MAPPING --------\n"
        "  1) Relationship between two people        (LCA on a DAG)\n"
        "  2) Everyone exactly k generations below   (BFS downwards)\n"
        "  3) All bloodlines back to origin ancestors (DFS upwards)\n"
        "  4) Degrees of separation between two people (BFS shortest path)\n"
        "  5) Export Graphviz .dot file\n"
        "  0) Exit\n";
}

static bool readInt(const std::string& prompt, int& out) {
    std::cout << prompt;
    std::string line;
    if (!std::getline(std::cin, line)) return false;
    std::istringstream ss(line);
    if (!(ss >> out)) {
        std::cout << "(please enter a whole number)\n";
        return true;   // stay in the menu
    }
    return true;
}

void runMenu(FamilyTree& t) {
    printMenu();
    while (true) {
        int choice = -1;
        if (!readInt("\nChoice> ", choice)) break;
        if (choice == 0) {
            std::cout << "Goodbye.\n";
            break;
        }

        if (choice == 1) {
            int a, b;
            if (!readInt("  Person A id: ", a)) break;
            if (!readInt("  Person B id: ", b)) break;
            if (!t.hasPerson(a) || !t.hasPerson(b)) {
                std::cout << "  ERROR: unknown id.\n";
                continue;
            }
            std::cout << "  " << t.relationshipText(a, b) << "\n";

        } else if (choice == 2) {
            int id, k;
            if (!readInt("  Person id: ", id)) break;
            if (!readInt("  How many generations below (k): ", k)) break;
            if (!t.hasPerson(id)) {
                std::cout << "  ERROR: unknown id.\n";
                continue;
            }
            std::vector<int> tier = t.generationsBelow(id, k);
            if (tier.empty()) {
                std::cout << "  Nobody exactly " << k
                          << " generation(s) below " << t.person(id).name
                          << ".\n";
            } else {
                std::cout << "  Exactly " << k << " generation(s) below "
                          << t.person(id).name << " (" << id << "):\n";
                for (std::size_t i = 0; i < tier.size(); ++i)
                    std::cout << "    - " << t.person(tier[i]).name << " ("
                              << tier[i] << ")\n";
            }

        } else if (choice == 3) {
            int id;
            if (!readInt("  Person id: ", id)) break;
            if (!t.hasPerson(id)) {
                std::cout << "  ERROR: unknown id.\n";
                continue;
            }
            std::vector<std::vector<int> > lines = t.lineagePaths(id);
            std::cout << "  " << lines.size() << " bloodline(s) from "
                      << t.person(id).name << " (" << id << "):\n";
            for (std::size_t i = 0; i < lines.size(); ++i)
                std::cout << "    " << (i + 1) << ". "
                          << t.lineageToString(lines[i]) << "\n";

        } else if (choice == 4) {
            int a, b;
            if (!readInt("  Person A id: ", a)) break;
            if (!readInt("  Person B id: ", b)) break;
            if (!t.hasPerson(a) || !t.hasPerson(b)) {
                std::cout << "  ERROR: unknown id.\n";
                continue;
            }
            std::vector<int> path;
            if (!t.shortestKinshipPath(a, b, path)) {
                std::cout << "  " << t.person(a).name << " and "
                          << t.person(b).name
                          << " are not connected (disconnected graph).\n";
            } else {
                std::cout << "  Degrees of separation: "
                          << (path.size() - 1) << "\n  chain: ";
                for (std::size_t i = 0; i < path.size(); ++i) {
                    if (i) std::cout << " - ";
                    std::cout << t.person(path[i]).name;
                }
                std::cout << "\n";
            }

        } else if (choice == 5) {
            std::cout << "  Writing family.dot ... "
                      << (t.exportDot("family.dot") ? "done." : "FAILED.")
                      << "\n  Render with: dot -Tpng family.dot -o family.png\n";

        } else {
            std::cout << "  Unknown choice.\n";
        }
    }
}

// ============================================================================
int main(int argc, char** argv) {
    FamilyTree t = buildSampleFamily();
    std::string mode = (argc > 1) ? argv[1] : "";

    if (mode == "--bench") {
        runBench();
        return 0;
    }
    if (mode == "--test") {
        runSelfTest(t);
        return gFailures == 0 ? 0 : 1;
    }
    if (mode == "--help" || mode == "-h") {
        std::cout << "Usage: family_tree [--test | --bench | --help]\n"
                     "  (no flag) : interactive menu\n"
                     "  --test    : validation demo + sample queries (self-checked)\n"
                     "  --bench   : timing experiment on random family DAGs\n";
        return 0;
    }

    std::cout << "Algorithmic Relationship Mapping via Graph Data Structures\n"
              << "Loaded " << t.size()
              << " people (4 generations). Type 0 to exit.\n";
    runMenu(t);
    return 0;
}

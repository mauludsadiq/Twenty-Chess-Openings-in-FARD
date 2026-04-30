# Twenty Chess Openings

A Fard language project that models 20 chess openings as an
information-theoretic channel: Opening -> Motif -> Endgame.

The repository quantifies uncertainty in the transition from opening
choice to resulting endgame type, using Shannon entropy measured over
a corpus of 325,977 master-level games.

## Repository Structure

data/
  openings.json           - 20 chess openings with ECO codes, SAN prefixes, tags
  endgames.json           - 20 theoretical endgame positions with material signatures
  motifs.json             - 20 transition motifs describing opening-to-endgame pathways
  Z_edges_fitted.json     - 22 directed edges with fitted probabilities (p_hat) and support counts
  Z_edges.json            - Raw/unfitted edges before probability estimation
  Z_edge_lengths.json     - Placeholder for edge path length statistics (not yet computed)

fard/
  apps/
    chess_query.fard      - Query an opening and return its probabilistic endgame distribution
    chess_stats.fard      - Compute entropies H(O), H(Z), H(E) and deltas
    chess_stats_report.fard - Full stats plus formatted endgame probability report
  packages/
    kernel/
      entropy.fard        - Shannon entropy utility function
  probes/
    probe_edge_fields.fard
    probe_fold_edges.fard
    probe_json.fard
    probe_runtime.fard
    probe_unique_openings.fard
  legacy/
    chess_stats_v0.fard   - Evolved through v0..v3 to the current stats program
    chess_stats_v1.fard
    chess_stats_v2.fard
    chess_stats_v3.fard

.gitignore                - Ignores /target and receipts/

## The 20 Openings

Alekhine's Defence, Scandinavian Defence, Dutch Defence, Philidor Defence,
Nimzowitsch Defence, Grünfeld Defense, Ruy Lopez Closed, Sicilian Moscow
Variation, Caro-Kann Advance, French Winawer, Closed Sicilian, Queen's Gambit
Accepted, Petroff Defense, 1.d4 e6 systems, Symmetrical English, Sicilian
Najdorf, Ruy Lopez Berlin, Slav Defense, Nimzo-Indian Rubinstein, King's
Indian Classical.

## The 20 Endgames

Lucena Position, Philidor Position, King Opposition, Rule of the Square,
Key Squares, Saavedra Position, Réti Maneuver, Wrong-Color Bishop, Lasker
Rook Defense, Vancura Position, Troitzky Line, Fortress, Shoulder Zugzwang,
Stamma Position, Kling-Horwitz Defense, Cochrane Defense, Two Knights vs Pawn,
Bishop vs Knight, Opposite-Color Bishops, Rook+Bishop vs Rook+Knight.

## The Graph

22 edges connect openings to endgames through motifs. Each edge has:
- p_hat: transition probability (fitted, some edges share probability with
  the same motif, e.g. Dutch splits 0.5/0.5 between Opposition and KeySquares)
- support_games: number of games in the corpus for that pathway
- motif_id: the named transition motif

Total support across all edges: 325,977 games.

## Information-Theoretic Results

H(O) = 4.002581 bits   (entropy of opening distribution)
H(Z) = 4.087176 bits   (entropy of joint opening-endgame-motif states)
H(E) = 4.023336 bits   (entropy of endgame distribution)

Delta H (emanation):  H(Z) - H(O) = 0.084595 bits
Delta H (collapse):   H(Z) - H(E) = 0.063840 bits

Adding motif and endgame detail to an opening adds about 0.085 bits of
surprise. Forgetting the motif and keeping only the endgame loses 0.064 bits.

## Usage

To query an opening's endgame distribution, edit chess_query.fard and
change the query_opening_id, then run with the Fard runtime.

## Data Sources

support_games counts are populated empirical counts used by the FARD kernel.
p_hat values are fitted conditional transition weights for opening-to-Z edges.
The legacy data_source label is stale metadata and does not control execution.

## Requirements

- Fard language runtime
- std libraries: fs, json, list, rec, math

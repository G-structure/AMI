### Background and Notations

Given a closed manifold triangle mesh M=(V,E)M=(V,E), a _seed_ vertex s∈Vs∈V, and a stitch width w∈Rw∈R, our goal is to generate human-readable instructions P(M,s,w)P(M,s,w)for crocheting MM from the point ss, with the given stitch width ww.

Crochet has a wide variety of stitches, and we focus here on the simple stitch used for Amigurumi, named _single crochet_ (sc). This is an approximately square stitch, thus covering MM with sc stitches is equivalent to constructing a quad re-mesh of MM, where each quad is a square, and all the edge lengths are constant. This is of course not possible unless the surface is developable, i.e. has zero Gaussian curvature. In practice, curved geometry is accommodated in crochet by introducing stitches which locally _increase_ (inc(x)) or _decrease_ (dec(x)) the amount of stitches by xx. Figure 3 (left) shows an example of the inc and dec stitches on a croched patch. Crochet instructions for Amigurumi typically include _rows_, where each row is a series of sc, inc, dec stitches. Figure 2 (middle) shows the instructions (pattern) for crocheting the sphere in Figure 2 (right).

### The Crochet Graph

A crochet stitch is composed of a _top_, a _base_ and a _stem_, where the inc, dec stitches have multiple stems, see Figure 3 (right). The top of one stitch is always the base of some stitch on the next row, and similarly, each stitch has a base on the previous row. Therefore, a natural abstraction of the stitch pattern is to consider the stitches and their interconnections as a _graph_.

Specifically, we define the _Crochet Graph_G=(S,R∪C)G=(S,R∪C), whose vertices SS are tops/bases of stitches, where a vertex (i,j)∈S(i,j)∈S is the base of the jj-th stitch in row ii, and the vertices in each row are consecutively ordered. The _column_ edges CC are stems of stitches, and the connectivity between the bases in each row is represented by the _row_ edges RR. We denote the total number of rows by NN. Figure 2 (left) shows the crochet graph corresponding to the croched sphere in Figure 2 (right).

A crochet graph is an intermediate representation between the input triangle mesh MM, and the output instructions PP. Our goal is to generate a graph such that it is (1) translatable to valid crochet instructions PP, and (2) when PP is crocheted and stuffed, the result resembles the input mesh. Note that there exist multiple instructions PP for the same graph GG, and within this space we aim for instructions which are _human-readable_.

We base our algorithm on the following observations.

Definition 2.1 ().: _A coupling (Gold and Sharir, 2018) $C=(c_{1},..,c_{k})betweentwosequencesbetweentwosequencesA=(p_{1},..,p_{n})andandB=(q_{1},..,q_{m})isanorderedsequenceofdistinctpairsofpointsfromisanorderedsequenceofdistinctpairsofpointsfromA\times B,suchthat,suchthatc_{1}=(p_{1},q_{1}),c_{k}=(p_{n},q_{m})$ and_

c_{r}=(p_{s},q_{t})\Rightarrow c_{r+1}\in{(p_{s+1},q_{t}),(p_{s},q_{t+1}),(p _{s+1},q_{t+1})},\quad\forall r<k. \tag{1}

Definition 2.2 ().: _Let $\mathcal{S}_{i},\mathcal{S}_{\mathrm{H}i},1\leq i<N,betheverticesoftwoconsecutiverowsof,betheverticesoftwoconsecutiverowsof\mathcal{G}=(\mathcal{S},\mathcal{R}\cup C),where,where\mathcal{S}_{i}=\big{(}(i,1),..,(i,n_{i})\big{)},where,where(i,j)\in\mathcal{S},and,andn_{i}isthenumberofverticesinrowisthenumberofverticesinrowi.Ifthereexistsacoupling.IfthereexistsacouplingCbetweenbetween\mathcal{S}_{i}andand\mathcal{S}_{\mathrm{H}i}suchthatforallsuchthatforallp_{s}\in\mathcal{S}_{i},q_{t}\in\mathcal{S}_{\mathrm{H}i+1}wehavethatwehavethat(p_{s},q_{t})\in Cifandonlyififandonlyif(p_{s},q_{t})\in C$, then the two rows are coupled._

Observation 2.3 ().: _If all the pairs of consecutive rows of GG are coupled, then there exist valid crochet instructions P(G)P(G) that use only the instructions sc, inc(x) and dec(x)._

Definition 2.4 ().: S→MS→M be an embedding of the vertices of GG on MM. An embedded edge of GG is a shortest geodesic between the embedding of two vertices of GG which share an edge, or between the embedding of the first and last vertices on the same row._

Definition 2.5 ().: _Let (p,q)∈M(p,q)∈M be two points whose geodesic distance is larger than some constant that depends on the stitch width ww, and let $p_{\mathcal{G},q}betheshortestgeodesicbetweenthem.Ifbetheshortestgeodesicbetweenthem.If\gamma_{p,q}$ intersects

Figure 3. (left) Crochet stitches inc,dec marked using yarn on a crocheted patch. (right) The anatomy of crochet stitches, marked are the _top_/_bottom_ of the stitch in red and the _stem_ in green.

Figure 2. (left) The _Crochet Graph_ of a sphere, with row edges RR in red, and column edges CC in blue. (middle) The corresponding instructions (pattern) for crocheting the sphere. (right) The crocheted sphere.

some embedded edge of an embedding XGXG​, for any two such points, then we say that XGXG​ covers MM._

Observation2.6.: _Let $X_{\mathcal{G}}:\mathcal{S}\to Mbeanembeddingoftheverticesofbeanembeddingoftheverticesof\mathcal{G}whichcoverswhichcoversM,and,andP(\mathcal{G})validcrochetinstructionsforvalidcrochetinstructionsfor\mathcal{G}.Ifalltheedgelengthsinducedby.IfalltheedgelengthsinducedbyX_{\mathcal{G}}areequaltoareequaltow,thenwhen,thenwhenP(\mathcal{G})willbecrochetedandstuffedtheresultwillbe"similar"towillbecrochetedandstuffedtheresultwillbe"similar"toM$._

We discuss Observation2.3 in section5.1, where we show how to translate the graph into valid instructions. The Observation2.6 is in fact true only for a subset of meshes, as we discuss in the next section.

### Crochetable Models

Curvature.Crocheting the patterns yields an empty flexible shell of fabric, which obtains its final shape by _stuffing_. Whether the stuffed model obtains the intended shape depends on how the model is stuffed (lightly, firmly, uniformly), as the yarn has some flexibility and will extend to accommodate if the model is overstuffed. We will assume that the model is stuffed enough to attain maximal volume, but not too much to cause the yarn to stretch and generate gaps. Thus, we expect the resulting stitch size to be similar to the edge lengths induced by the embedding of the crochet graph. If for a given graph GG its embedding in 3D with edge lengths ww is _unique_, we expect the crocheted and stuffed shape to be similar to the input surface.

Importantly, unless the shape is convex, the edge lengths alone (i.e., the _metric_) do not contain enough information to uniquely determine the shape of a non-stuffed model. For example, a surface that has a "crater" leads to edge lengths which can be realized either as a crater or as a "hill". However, if we add the maximal volume assumption, only the hill is possible. This in fact implies that surfaces which have "craters", or more formally, regions of negative mean curvature with _positive_ Gaussian curvature, cannot be realized by crocheting and stuffing alone. This is similar to the observation made by Konakovich et al. (2018), that surfaces which have negative mean curvature cannot be realized by maximizing the volume with a given conformal metric (i.e., only isotropic scaling is allowed relative to the flat configuration). We handle this case similarly, by preprocessing the model so that it does not contain "craters" (Section7.1.1).

Furthermore, in our case, since we allow anisotropic scaling, negative mean curvature with _negative_ Gaussian curvature is in fact possible, but requires a modified sampling rate, which we discuss in Section7.1.2. To conclude, in terms of geometric obstructions to crochetibility, the four possible curvature situations are summarized in Table1. Note that, as with any sampling-dependent method, if the sampling rate (namely, the number of rows NN) is too small compared to the feature size of the model, the crocheted output will lose some of the geometric detail.

Branching.Observation2.6 requires that the crochet graph embedding XGXG​ covers the input surface MM. Because of the special structure of this graph, this induces additional constraints on the possible geometries. Intuitively, models that branch (see Figure7) cannot be covered in this way. Mathematically, this means that the geodesic distance function on MM from the embedding of the seed vertex ss cannot have saddles. This is solved by segmenting the shape, and crocheting the segments in an iterative manner. We explain this in detail in Section7.2.

We first explain the generation of the crochet pattern for a simple non-branching model with positive mean curvature, and then discuss how we handle negative mean curvature and branching.

## 3. Overview

Given a 3D mesh MM, a seed point ss and a stitch width ww, we first compute a crochet graph GG and its embedding XGXG​ such that they adhere to Observations2.3 and 2.6 (Section4). Then we compute the crochet pattern P(G)P(G) (Section5).

To generate GG and XGXG​, we first compute the embedding of the vertices SS on MM (Section4.1), and then derive from that the connectivity of GG, i.e. the row edges RR and column edges CC (Section4.2). To compute the pattern P(G)P(G), we first translate the graph into a program using standard code synthesis tools (Section5.1), and then apply loop unrolling to make the pattern human-readable (Section5.2). See Algorithm1.

## 4. Mesh to crochet graph

### Geometry

Observation2.3 implies that the vertices SS should be grouped into ordered rows, where in each row the vertices have a well defined order. We address this requirement by computing two non-negative monotonically increasing, constant speed functions f,g:M→Rf,g:M→R which define the row-order and column-order of every point on MM. Furthermore, Observation2.6 implies that the distance between embedded rows, and between embedded vertices in the same row should be ww. We address this by sampling f,gf,g appropriately.

Row order ff.Our models are closed (so they can be stuffed), and therefore the first and last rows in the graph GG contain a single

```
Input: A triangle mesh $M$, seed $s$, stitch width $w$ Output: Embedded crochet graph $\mathcal{G}=(\mathcal{S},\mathcal{R}\cup\mathcal{C}),X_{\mathcal{G}}$, crochet pattern $P(\mathcal{G})$ Mesh to Graph :// Section4  Geometry $\mathcal{S},X_{\mathcal{G}}$ :// Section4.1  Connectivity $\mathcal{R},\mathcal{C}$ :// Section4.2 Graph to pattern :// Section5  Graph to program :// Section5.1  Program to pattern $P(\mathcal{G})$ :// Section5.2
```

\begin{table} \begin{tabular}{c|c|c} Mean & & \ Curvature & & \ Gaussian & Positive & Negative \ Curvature & & \ \hline Positive & Crochetable & Preprocessing \ Negative & Crochetable & Sampling modification \ \end{tabular} \end{table} Table 1. Curvature obstructions to crochetability, see the text for details.

vertex. The first row contains only the seed ss, and its function value is f(s)=0f(s)=0. We take f(v),v∈Vf(v),v∈V to be f(v)=d(v,s)f(v)=d(v,s), where dd is the geodesic distance. Thus, the isolines of ff are rows, and two points p,q∈Mp,q∈M are on the same row if f(p)=f(q)f(p)=f(q). If ff has more than one maximum, then we need to handle branching (see Section 7.2). Otherwise, the vertex that attains the maximum of ff, denoted as fMfM​ will be the single vertex on the last row.

Column order ggWe first cut MM along a geodesic from ss to fMfM​, so that our model and the graph that we compute have the same topology, and denote the cut model by MCMC​. The requirements are that _within each row_ the vertices of GG have a well defined order. A row is an isoline of ff, and therefore the rate of change along the isoline is given by the directional derivative of gg in the direction of the tangent to the isoline. Specifically, the tangent to the isoline of ff at a point p∈Mp∈M is given by J∇fJ∇f, where ffis the rotation by π/2π/2 in the tangent plane of pp. Thus to find gg, we solve an optimization problem whose objective is to minimize ∫MC∣⟨∇f,∇g⟩−1∣2∫MC​​∣⟨∇f,∇g⟩−1∣2, s.t. g(B)=0g(B)=0. Here, B⊂VB⊂V is the longest connected path of boundary vertices of MCMC​ along which ff is strictly monotone.

SamplingThe functions f,gf,g define a parameterization of MCMC​ to the plane. We conjecture that this parameterization is bijective (as it was in all of our experiments), but leave the proof to future work. The parameterization may have a large metric distortion, however, if f(p)=f(q)=f0f(p)=f(q)=f0​ for some two points p,q∈Mp,q∈M, then ∣g(p)−g(q)∣∣g(p)−g(q)∣ is equal to the length of the isoline of f0f0​ between pp and qq. Therefore, we uniformly sample f,gf,g on a 2D2D grid of width ww, yielding the vertices of SS with indices (f/w,g/w)(f/w,g/w). Pushing forward the sampled points to the mesh MCMC​ yields the embedding of SS on MCMC​ (and therefore MM), namely XGXG​.

### Connectivity

Row edgesRREach two consecutive vertices of SS on the same row are connected by a row edge. Namely, $\mathcal{R}=\bigcup_{j=1}^{N}\mathcal{R}_{i},and,and\mathcal{R}_{i}=\left{(i,j),(i,j+1)\right}|\ j\in{1,\neg,n_{i}-1}\right}.Here.Heren_{i}=\left|\mathcal{S}_{i}\right|=\left|{(i,j)\in\mathcal{S}}\right|,namelythenumberofverticesinthe,namelythenumberofverticesinthei$-th row.

Let x,y∈Sx,y∈S be two consecutive vertices on the ii-th row. Then we have that f(x)=f(y)=f0f(x)=f(y)=f0​ and ∣g(x)−g(y)∣=w∣g(x)−g(y)∣=w. Therefore, dγ(f(f))(XG(x),XG(y))=wdγ(f(f))​(XG​(x),XG​(y))=w, where γ(f0)γ(f0​) is the isoline of f0f0​ on MM, and dγ(f0)dγ(f0​)​ is the distance along the isoline. Hence, the Euclidean distance between the embedded vertices ∣∣XG(x)−XG(y)∣∣≤w∣∣XG​(x)−XG​(y)∣∣≤w, and the distance tends to ww for a "small enough" stitch size. Here, "small enough", means on the order of the square root of the radius of curvature of γ(f0)γ(f0​), which is given in terms of the normal curvature in direction J∇fJ∇f.

Column edgesCCFirst, Observation 2.3 requires that all pairs of consecutive rows are coupled. Let CiCi​ be the coupling corresponding to rows $\mathcal{S}_{i},\mathcal{S}_{i+1},andlet,andlet(p_{\mathrm{s}},q_{t})\in C_{i}.Since.Sincep_{\mathrm{s}}andandq_{t}areonconsecutiverows,andthereforeembeddedonisolinesofareonconsecutiverows,andthereforeembeddedonisolinesoffwhichdifferbywhichdifferbyw,theminimaldistance,theminimaldistance||X_{\mathcal{G}}(p_{\mathrm{s}})-X_{\mathcal{G}}(q_{t})||isclosetoisclosetow$. Therefore, if among all couplings we seek the minimizer of:

\min_{C_{i}coupling}\sum_{(p_{\mathrm{s}},q_{t})\in C_{i}}||X_{\mathcal{G}}(p_{ \mathrm{s}})-X_{\mathcal{G}}(q_{t})||, \tag{2}

then the length of the column edges will be close to ww.

A minimal coupling between every pair of consecutive rows is found by Dynamic Time Warping (DTW) (Gold and Sharir, 2018; Sakoe and Chiba, 1978).

## 5. Crochet graph to instructions

### Graph to Program

In order to turn the crochet graph into instructions, we rely on the following observation: crochet instructions (patterns) constitute an _instruction set_, and as such, crocheting is an _execution_ and the finished object is an _execution result_. Moreover, because of the nature of a crocheted object, it is not only a result but a step-by-step _execution trace_.

Therefore, given a crochet object, or in this case its graph representation GG, deriving instructions constitutes a form of _execution reconstruction_(Zuo et al., 2021), a method for reconstructing the set of instructions that lead to an execution trace. While reconstructing an execution trace generally requires searching an exponential space of possible instruction sequences, the crochet instruction set

Figure 4. (a) The isolines of ff, with the seed (red), the maxima of ff (blue), and the cut (black) marked. (b) The f,gf,g parameterization, and the sampled grid SS, (c) The pushed forward points XGXG​. (d) The output crochet graph GG, with red row edges RR and blue column edges CC.

[MISSING_PAGE_FAIL:6]

et al., 2012) localized to these areas until the mean curvature is positive everywhere.

#### 7.1.2. Negative Gaussian curvature

The sampling rate of the isolines of ff is determined by the directional derivative of gg w.r.t. the tangent to the isoline, namely by ⟨∇g,J∇f⟩⟨∇g,J∇f⟩. If the curvature in the direction of the isoline kfvfkf​vf​ is large compared to the stitch width ww, a uniform sampling rate is inadequate, and does not result in a similar geometry. We therefore adjust the sampling rate in these regions by setting ⟨∇g,J∇f⟩=h(kfvf)⟨∇g,J∇f⟩=h(kf​vf​). We take h(x)=tanh⁡(−x/α)/2+1h(x)=tanh(−x/α)/2+1, to avoid degenerate sampling rates, with α=10α=10. Figure 6 shows an example of such a model. In the central region, the model has negative mean and Gaussian curvature (a,b). Using a uniform sampling rate does not fully reconstruct the negative curvature (d,f), whereas a curvature adapted sampling does (c,e).

### Branching

If for a given seed ss, the geodesic function f(x)=d(s,x)f(x)=d(s,x) has multiply-connected isolines, the graph GG cannot cover the model. In these cases, the model is automatically decomposed into multiple segments, each of which can be crocheted using the approach described in Sections 4 and 5. The segments are attached by a method called "join-as-you-go" (Bennett, 2020), meaning each segment is crocheted onto the last row of the previous segment, and therefore no additional sewing is required. Furthermore, the segment boundaries are not visible in the crocheted object. The more common method, in contrast, involves sewing together closed segments, which requires accuracy to achieve geometric details such as symmetry. In this section we describe the modifications required to accommodate such objects.

#### 7.2.1. Mesh to Graph

Given a seed ss let f(x)=d(s,x)f(x)=d(s,x). Let Π=(σ1..,σm)Π=(σ1​..,σm​) be the saddle points of ff, sorted by fi=f(σi)fi​=f(σi​). Namely f1≤f2≤...≤fmf1​≤f2​≤...≤fm​. For each σiσi​ in order, we compute the isoline of fifi​, denoted by yiyi​, and slice a new segment for each connected component of yiyi​. Meaning, the segments are obtained by slicing along the isolines of the saddle points, ordered by increasing geodesic distance to the seed. Figure 7(a) shows the isolines of ff for the seed

Figure 6. **A model with negative mean curvature (a) and negative Gaussian curvature (b) in the same region. The corresponding embeddings YGYG​ (c,d) and knitted objects (e,f), with curvature-adapted (c,e) and uniform (d,f) sampling rate of gg. Note that the uniform sampling rate does not lead to a model that is similar to the input, whereas curvature-adapted sampling yields a better result.**

Figure 7. **(a) The isolines of ff for the marked red seed point, and the saddles (cyan) and maxima (blue) of ff. (b) The resulting segments. (c) The crochet graph G,XGG,XG​. (d) The embedding YGYG​ of the crochet graph. (e) The final knitted model. Different segments were crocheted in different colors for better visualization.**

point ss marked in red, as well as the saddles (cyan), and maxima (blue).

In addition to the segmentation, we generate a directed graph GσGσ​, whose vertices are the segments, and where an edge (s,t)(s,t) exists if the segments MsMs​ and MtMt​ share a boundary, and ff values on MsMs​ are smaller than ff values on MtMt​. The crocheting order of the segments is determined by a topological sort of GσGσ​. Figure 7(b) shows the resulting segments. Very thin segments might not be sampled (marked in black in Figure 7(b) and other figures), and are skipped and not crocheted.

_Geometry._ Any resulting segment MfMf​ is either a half sphere or a cylinder, and thus can be covered by a crochet graph GfGf​. While ff is computed for the whole model before segmentation, gg is computed for each segment separately. The cut for gg is made from the maxima of ff in the segment to the closest point on the segment's boundary. If ffattains its maximum on one of the boundaries of the segment (there are at most two boundaries), then cut is computed to the closest point on the other boundary.

Figure 8 shows an example of the location of the cut, where we show the front and the back of the model. We show the location of the cut in brown on the segmented model (a), as well as the location of the cut in the crocheted model (b). To mark the location of the cut during crocheting a piece of yarn was used to mark the beginning/end of the row.

_Connectivity._ For every two segments Ms,MfMs​,Mf​ which share an edge (s,t)(s,t) in GσGσ​, we add an additional condition that the last row of MsMs​ is coupled to the first row of MfMf​. Figure 7(c) shows the crochet graph for all the segments of the Homer model. Finally, (e) shows the crocheted model, where each segment was crocheted with a different color for better visualization.

#### 7.2.2. Graph to instructions.

The same simulation of the crochet operations is applied to the first row of a new segment, but the sequence of stitches that is used as its previous row is no longer a full row. Instead, the last rows of all attached segments are arranged and filtered to include only vertices that have a connecting edge to the new segment's row, constituting a "joint" previous row. The transducer then takes stock of when its consumption of the previous row skips stitches, splits segments, skips segments, or spans multiple parent segments, and includes this information in the row instructions.

## 8. Limitations

Our method only handles closed surfaces, since we aim for Amigurumi models, which are stuffed. The segmentation approach may generate thin segments, which are harder to crochet. While we filter out very thin segments, we believe that small modifications to the singularity locations can yield a better segmentation without considerably affecting the shape. Our approach does not take the symmetry of the model into account, and thus discretization errors for low resolution patterns may lead to non-symmetric crocheted models. Our current setup also has limited design freedom. While using a single seed is simple, it does not give the user control of the knitting direction throughout the shape.

## 9. Results

### Implementation Details

We implemented our algorithm in Matlab and C++. We use Heat Geodesics (Crane et al., 2017) for computing geodesic distances, where we took the recommended time parameter tt, namely the average edge length squared. In some cases, the mesh has too many geometric details, leading to neighboring saddles/ extrema of the distance function. In this case, we take advantage of the tunability of heat geodesics, and repeated multiply tt by an increasing power of 22 until there are no more neighboring saddles/extrema. For computing geodesic paths we use (Sharp and Crane, 2020). Table 2 provides the statistics for all the models. Algorithm running times were measured on a desktop machine with an Intel Core i7.

Figure 8. (a) The cut location (brown) on a segmented model. (b) The crocheted model with the cut marked using a piece of yarn.

### Gallery

Figures 1, 7, 6 demonstrate our results. Figure 10 shows additional models, where we show the (a) segmented shape and seed, (b) the crochet graph GG and its embedding XGXG​ on the input mesh, (c) the embedding YGYG​ generated from the edge lengths, and (d) the final crocheted object. The sixth row of Figure 10 shows an example of a pretzel model, which cannot be scheduled for machine-knitting as discussed at AutoKnit (Narayanan et al., 2018) (see Figure 25 there). Our method, on the other hand, generates crochetable instructions. The last two rows of Figure 10 show the results for the same model, seed point, yarn type and hook but using different stitch width ww. Note that while manually adapting the pattern to different sizes is a difficult task, our algorithm preforms it automatically. Note that the shapes are similar to the input, and the segment boundaries are not visible in the output crocheted models. Therefore, we can achieve varied geometries without visibly segmenting the shape, leading to more visually pleasing results than the approach by Igarashi et al. (2008a) (see Figure 11 there). Furthermore, our results have a closer resemblance to the input, compared to the method by Guo et al. (2020) (see Figure 12 there).

### Creases

Since stuffed items tend to be smooth, there exist crochet shaping techniques that allow the generation of creases. Specifically, instead of inserting the hook under both loops of the stitch, the hook is inserted only in the front loop (denoted _Front Loop Only_ FLO), or only in the back loop (denoted _Back Loop Only_ BL0). The BL0 (resp. FLO) stitch allows creases which are positively (resp. negatively) curved with respect to the columns direction (i.e. orthogonal to the knitting direction).

We define _crease vertices_ as vertices that have large maximum absolute curvature, and their maximum absolute curvature direction is orthogonal to the knitting direction. For any two consecutive crease vertices on the same row, we mark all the stitches based on these vertices as BL0 or FLO, depending on the sign of the curvature. We allow the user to choose whether to enable this option or not. Figure 9 shows an example of a model where this shaping technique was used. The BL0 (resp. FLO) edges are marked in green (resp. cyan) in (b). Note the corresponding sharp crease in the knitted object (c). We note that creases have also been used in knitting. For example, (Wu et al., 2019) allow for creases, but they use a different technique.

## 10. Conclusions and Future Work

We presented a novel automatic approach for generating crochet knitting instructions for Amigurumi from an input triangle mesh. Given a single seed point and a stitch size, we generate human-readable instructions that use only simple crochet stitches (sc, inc, dec). Our method is applicable to a variety of geometries, and leads to crocheted models which are visually similar to the input shape. In the future we plan to incorporate our method within an interactive framework that allows the user to move the seed point, change the yarn and the gauge and see the expected shape of the output. Furthermore, we plan to add colors and texture, as well as support for additional types of stitches.

With the wide popularity of Amigurumi, and crochet in general, in recent years, we believe that tools that allow novice users, as well as pattern designers, to generate crochet instructions from 3D models would be quickly adopted and built upon by the crocheters and makers communities. Our approach provides an important stepping stone in this direction, and we expect that it will sow the seeds for further research in computational crochet.

## Acknowledgments

M. Edelstein acknowledges funding from the Jacobs Qualcomm Excellence Scholarship and the Zeff, Fine and Daniel Scholarship. M. Ben Chen acknowledges the support of the Israel Science Foundation (grant No. 1073/21), and the European Research Council (ERC starting grant no. 714776 OPREP). We also thank SHREC'07 and Josh Holianty for providing models.

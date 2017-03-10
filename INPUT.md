STATIS est une méthode exploratoire d’analyse des données utilisée pour les données quantitatives. [Escouffier et L’Hermier des Plantes, 1976]. Elle réalise l’exploration simultanée de plusieurs tableaux de données recueillis sur les mêmes individus. Le nombre de variables peut différer d’un tableau à un autre. On parle communément de T études $(X_{t}, M_{t}, D)_{t=1,...,T}$.\newline
  Il est à noter qu'une étude est le triplet (X, M, D), où:
    \begin{itemize}
    \item X est le tableau de données centré croisant individus et variables;
    \item M est la métrique  permettant le calcul de distances entre individus, de taille $n*n$;
    \item D est la métrique des poids, permettant le calcul de distances entre variables, de taille $p*p$.
    \end{itemize}

  Le but de STATIS est la recherche d’une structure commune aux études, appelée infrastructure. Pour la mettre en place, il faut passer par les quatre étapes successives  suivantes :  
\section{Étude de l’interstructure}
  On y compare globalement les T tableaux de données. L'objet $W_{t}$ est représentatif de l'étude (t). Il est exprimé par la formule suivante :
  \begin{equation}
  W_{t} = X_{t}M_{t}X'_{t} 
  \label{Wt}
  \end{equation}
  Le produit scalaire de Hilbert-Schmidt est utilisé pour définir une distance entre objets $W_{t}$.
  \begin{equation}
  \langle W_{t}|W_{t'}\rangle_{HS}=Tr(DW_{t}DW_{t'})\label{produitHS}
  \end{equation}
  A partir de ce produit est définie la norme de l’objet $W_{t}$ et la distance entre les objets $W_{t}$ et $W_{t'}$:
  \begin{equation}
  ||W_{t}||^{2} = \big\langle X_{t}|X_{t}\rangle_{HS} =  \sum_{l=1}^{n} (\lambda_{l}^{(t)})^{2}
  \end{equation}
  \begin{equation}
  d_{HS}(W_{t},W_{t'})=||W_{t}-W_{t'}||_{HS}=\sqrt[]{||W_{t}||_{HS}^{2}+||W_{t'}||_{HS}^{2}-2\langle X_{t}|X_{t'}\rangle_{HS}}
  \end{equation}

  Dans le cas où les objets ont des normes trop différentes, il est primordial de les normaliser en utilisant l'objet normé $W_{t}/||W_{t}||_{HS}$.
  \newline
  Le produit scalaire permet de déterminer la matrice des produits scalaires entre objets $W_{t}$, appelée S. Sa diagonalisation permet d'obtenir une image euclidienne des T études.
  
  Soit RV, le coefficient d'association qui est une corrélation vectorielle entre études. Il est utile pour l'interprétation de l'interstructure.

  \begin{equation}
  RV(t,t') = \Bigg \langle \frac{W_{t}}{||W_{t}||_{HS}} \Bigg| \frac{W_{t'}}{||W_{t'}||_{HS}} \Bigg \rangle = \frac{S_{tt'}}{\sqrt{S_{tt} S_{t't'}}}
  \end{equation}

  La valeur de RV varie entre 0 et 1. Plus les deux matrices $W_t$ et $W_t'$  sont similaires et plus sa valeur s'approche de 1.\newline
   
  L'image euclidienne des T études est obtenue, après avoir affecté un poids $\pi_{t}$ à chaque étude. Elle sert à visualiser l'interstructure. Soit $\Delta$ une matrice diagonale qui contient tous les poids $\pi_{t}$. \newline 
  Après avoir effectué une ACP de la matrice S, on obtient les points $A_{1},..., A_{t}$, qui sont associés aux études $W_{1},..., W_{t}$ dans l'image euclidienne obtenue. Dans ce registre:
  \begin{itemize}
  \item $\tau_{i}$ est la $i^{\grave{e}me}$ valeur propre de la matrice $S\Delta$ associée à $\gamma_{i}$, son $i^{\grave{e}me}$ vecteur propre.
  \item Les coordonnées des $A_{t}$ sur l'axe i sont les composantes du vecteur $\sqrt[]{\tau_{i}}\gamma_{i}$
  \end{itemize}
  La distance entre deux points $A_{t}$  et  $A_{t'}$ est la meilleure approximation possible de la distance de HS entre les objets $W_{t}$ et  $ W_{t'}$. Il est à noter que dans ce cadre, $RV(t,t')$ représente le cosinus de l'angle entre les vecteurs $OA_{t}$ et $OA_{t'}$.

  Comme dans l'ACP, si l'on ne veut pas que des tableaux interviennent dans la constitution de l'image euclidienne, on les ajoute comme tableaux supplémentaires.


\section{Recherche d'un compromis}
  Un compromis est une moyenne pondérée entre les objets $W_{t}$.
  \begin{equation}
  W=\sum_{t=1}^{T}\alpha_{t}W_{t}
  \end{equation}
  Les coefficients $\alpha_{t}$ sont déterminés en respectant les deux critères suivants : 

  \begin{itemize}
  \item Le compromis W est l’objet le plus corrélé en moyenne avec les objets $W_{t}$.
  \item L'objet W est de même nature que les objets $W_{t}$ (Sa norme étant la moyenne pondérée des normes des objets).
  \end{itemize}

  En considérant $\gamma_{1}$ le vecteur propre de S associé à la plus grande valeur propre $\lambda_{1}$, les coefficients $\alpha_{t}$ sont déterminés comme suit :
  \begin{equation}
  \alpha_{t}=\frac{1}{\sqrt[]{\lambda_{1}}} \Big ( \sum_{\tau=1}^{T}\pi_{\tau} \sqrt[]{S_{\tau\tau}} \Big ) \pi_{\tau}\gamma_{1}^{(t)}
    \end{equation}
 
 \begin{equation}
  W = \sum_{t=1}^{T} \Bigg[\Bigg(\frac{1}{\sqrt[]{\lambda_{1}}} \bigg(\sum_{\tau=1}^{T}\pi_{\tau} \sqrt[]{S_{\tau\tau}} \bigg) \pi_{\tau}\gamma_{1}^{(t)} \Bigg) W_{t} \Bigg]
  \end{equation}
  Dans l’image euclidienne des études, le compromis est situé sur le premier axe.

\section{Étude de l’intrastructure}
  L’ACP de la matrice WD fournit l’image euclidienne compromis.
  Soient $\mu_{1},..., \mu_{n}$ les valeurs propres de la matrice WD associés aux vecteurs propres $\varepsilon_{1},..., \varepsilon_{n}$. L’image euclidienne compromis des individus est composée des points $B_{1},...., B_{n}$; leurs coordonnées sur l’axe k sont les composantes du vecteur :
  \begin{equation}
  \sqrt[]{\mu_{k}} \varepsilon_{k} = \frac{1}{\sqrt[]{\mu_{k}}}WD\varepsilon_{k}
  \end{equation}
  
  L’étude de l’interstructure montre l’existence d’une structure des individus commune aux études, la représentation de l’image euclidienne compromis permet sa description.
  La distance entre deux points $B_{i}$ et $B_{j}$ de l’image euclidienne compromis représente la distance compromis entre les individus i et j. Elle est aussi la distance moyenne entre i et j sur la période étudiée.\newline
  
  Il est possible d’interpréter la position des individus sur un axe quelconque, noté k. Pour cela, on calcule les corrélations de la composante principale du compromis correspondant à cet axe avec les variables de chaque étude.
  \begin{equation}
  \big \langle \varepsilon_{k}, (x^{j})^{(t)} \big \rangle_{D} = (x^{j})^{(t)} \prime D\varepsilon_{k}
  \end{equation}
  Ces corrélations peuvent être résumées sur un graphique. L’étude de ce dernier est utile pour expliquer les positions compromis des individus dans leur image euclidienne.\newline
  Il est à noter que l'image euclidienne compromis est identique à celle trouvée avec une ACP sur une juxtaposition des tableaux initiaux multipliés par $\sqrt[]{\alpha_{t}}$.

\section{Représentation des trajectoires des individus}
  Les trajectoires décrivent les écarts des objets entre eux et avec le compromis au niveau individuel. Leur représentation se fait dans l’image euclidienne compromis en même temps que les T nuages d’individus débouchant sur une représentation de nT points.
  Pour obtenir la trajectoire, on place les différentes positions d’un individu tel qu’il est décrit par chaque étude (t). On se base sur les coordonnées des points compromis $B_{1},..., B_{n}$ sur l’axe k.
  En considérant chaque étude comme étant placée en supplémentaire, les coordonnées des points $B^{t}_{1},..., B^{t}_{n}$ sont pour $t=1...T$ :
  
  \begin{equation}
  \frac{1}{\sqrt[]{\mu_{k}}}W_{t}D\varepsilon_{k}
  \end{equation}
  
  Ces points n’ont pas influencé la construction de l’image euclidienne compromis, mais ils ont pu y être représentés.
  Les trajectoires décèlent les individus i responsables des écarts entre les études (t) et (t’). Elles s’interprètent par rapport à l’évolution moyenne. Et il y a deux grandes familles de formes de trajectoires:
  \begin{itemize}
  \item Peu étendue tournant autour de sa position compromis ce qui correspond à un individu suivant une évolution moyenne. 
  \item Ou de grande amplitude reflétant un changement de structure de l’individu au cours des années.
   \end{itemize}
  Si un individu est absent de certaines études, on le place dans l’image euclidienne compromis en le traitant comme un individu supplémentaire et en calculant sa trajectoire.
  Les individus placés en supplémentaire n’interviennent pas dans la détermination de l’infrastructure ni du compromis, ils n’ont pas de position compromis.
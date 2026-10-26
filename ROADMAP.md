# treeshapesim: estado e próximos passos

Atualizado em 2026-10-25, depois do merge do PR \#3 (`main` em
`4a65352`, versão de desenvolvimento 0.1.0.9000). A CI do `main` (R CMD
check em macOS, Windows e Ubuntu, e pkgdown) passou.

## Feito

| Item | O que entrou | Onde |
|----|----|----|
| Validação e robustez | `n`, `beta` e `c` validados (`c = Inf` dava ramos `NaN`); topologia construída sem recursão (n = 3.000 com β = −1,99 deixou de falhar), com as mesmas árvores para a mesma semente | PR \#2 |
| R CMD check na CI | workflow `R-CMD-check.yaml` (r-lib/actions v2) e badge no README | PR \#2 |
| 1\. Simular para valores-alvo | [`calibrate_shape()`](https://diogoprov.github.io/treeshapesim/reference/calibrate_shape.md) (curvas γ(c) e Colless(β), monótonas por regressão isotônica), [`shape_params()`](https://diogoprov.github.io/treeshapesim/reference/shape_params.md) (inversão) e [`sim_tree_target()`](https://diogoprov.github.io/treeshapesim/reference/sim_tree_target.md) (rejeição separada de topologias e intervalos, equivalente à rejeição de árvores inteiras) | PR \#3 |
| 2\. Grades de formato | [`shape_grid()`](https://diogoprov.github.io/treeshapesim/reference/shape_grid.md): cruza níveis de γ e Colless e dá β e `c` de cada célula | PR \#3 |
| 3\. Várias árvores por chamada | `nsim` em [`sim_tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/sim_tree_shape.md), que devolve `multiPhylo` | PR \#3 |
| 4\. Profundidade da raiz | `root_age` em [`sim_tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/sim_tree_shape.md) e [`sim_tree_target()`](https://diogoprov.github.io/treeshapesim/reference/sim_tree_target.md) | PR \#3 |
| 5\. Desempenho | **Parcial.** [`assign_ages()`](https://diogoprov.github.io/treeshapesim/reference/assign_ages.md) ficou sem [`ape::prop.part()`](https://rdrr.io/pkg/ape/man/boot.phylo.html) (ganho pequeno, árvores idênticas). O gargalo é [`rtopology_beta()`](https://diogoprov.github.io/treeshapesim/reference/rtopology_beta.md), O(n²) em árvores desbalanceadas; pré-calcular os `lgamma` mudava algumas árvores por arredondamento e foi revertido | PR \#3 |
| 6\. Vignette | [`vignette("treeshapesim")`](https://diogoprov.github.io/treeshapesim/articles/treeshapesim.md): separação dos processos, calibração, alvos e grade | PR \#3 |
| 7\. Citação | `inst/CITATION` (ainda sem DOI) | PR \#3 |
| 11\. [`tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/tree_shape.md) | γ e Colless normalizado (Yule) em R puro, iguais a [`ape::gammaStat()`](https://rdrr.io/pkg/ape/man/gammaStat.html) e `treestats::colless(normalization = "yule")`; o treestats continua em Suggests | PR \#3 |
| 13\. Validar [`assign_ages()`](https://diogoprov.github.io/treeshapesim/reference/assign_ages.md) | intervalos positivos e finitos; árvore binária e enraizada | PR \#3 |
| 14\. Documentação precisa | a independência de γ e Colless em relação a `c` e β agora é descrita como exata (para um dado n) | PR \#3 |

## A fazer

### Alta prioridade

- **DOI do Zenodo.** Ligar o repositório ao Zenodo e publicar uma
  release (por exemplo, 0.2.0), depois acrescentar o DOI ao
  `inst/CITATION` e ao README. Fecha de vez o item 7.
- **Versão 0.2.0.** Trocar 0.1.0.9000 por 0.2.0 no DESCRIPTION e no NEWS
  ao publicar a release, já que entraram funções novas.
- **Desempenho da topologia (resto do item 5).** Só vale se alguém
  precisar de árvores muito grandes e desbalanceadas. Para não mudar as
  árvores existentes, o caminho seria um gerador novo e opcional (por
  exemplo, em C++ ou com outro algoritmo de amostragem), mantendo o
  atual como padrão.

### Média prioridade

- **8. Amostragem incompleta e extinção.** Um argumento para a fração
  amostrada (podar pontas ao acaso) ou um processo de nascimento e morte
  com poda. O γ é sensível à amostragem incompleta.
- **9. Taxa de especiação definida pelo usuário.** `rate = function(k)`,
  mantendo k^c como padrão.
- **10. Guardar os parâmetros na árvore.** Gravar `beta`, `c` e
  `root_age` como atributos do `phylo` em
  [`sim_tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/sim_tree_shape.md);
  [`sim_tree_target()`](https://diogoprov.github.io/treeshapesim/reference/sim_tree_target.md)
  já guarda os seus.
- **12. Aviso para `c` extremo.** Valores finitos muito grandes (por
  exemplo,
  400. geram ramos de comprimento zero, sem aviso.
- **Calibração em cache.** Guardar curvas pré-calculadas para tamanhos
  comuns (por exemplo, 50, 100 e 200 pontas), para evitar a calibração a
  cada sessão. Exige uma rodada de simulação só para isso.

### Baixa prioridade

- **15. Outros modelos de topologia**, como o modelo alpha de Ford. As
  referências ainda não foram conferidas.
- **16. CRAN.** Faz mais sentido depois do aceite do artigo; exige
  `cran-comments.md`, as referências no campo Description no formato
  `<doi:...>` e checks no win-builder e no R-hub.
- **17. Cobertura de testes** (`covr` e badge).
- **18. Revisão no rOpenSci.**
- **19. Partes em C++**, se o desempenho continuar sendo um problema.

## Ligação com o artigo

- O `pcps-sim` usa `R/trees.R`, que gera as mesmas árvores que
  [`sim_tree_shape()`](https://diogoprov.github.io/treeshapesim/reference/sim_tree_shape.md)
  com a mesma semente; o teste com valores congelados protege isso.
- Depois do aceite: citar o treeshapesim nos Métodos ou no Apêndice S1,
  e só então ligar o pacote ao artigo no README e no site (anonimato da
  revisão).

# treeshapesim: estado e próximos passos

Atualizado em 2026-11-05, depois dos itens 8, 9, 10, 12 e 17 e da
correção de `c` para amostragem incompleta. Última versão publicada:
0.2.1; versão de desenvolvimento 0.2.1.9000.

## Feito

| Item | O que entrou | Onde |
|----|----|----|
| Validação e robustez | `n`, `beta` e `c` validados (`c = Inf` dava ramos `NaN`); topologia construída sem recursão (n = 3.000 com β = −1,99 deixou de falhar), com as mesmas árvores para a mesma semente | PR \#2 |
| R CMD check na CI | workflow `R-CMD-check.yaml` (r-lib/actions v2) e badge no README | PR \#2 |
| 1\. Simular para valores-alvo | [`calibrate_shape()`](https://provetelab.org/treeshapesim/reference/calibrate_shape.md) (curvas γ(c) e Colless(β), monótonas por regressão isotônica), [`shape_params()`](https://provetelab.org/treeshapesim/reference/shape_params.md) (inversão) e [`sim_tree_target()`](https://provetelab.org/treeshapesim/reference/sim_tree_target.md) (rejeição separada de topologias e intervalos, equivalente à rejeição de árvores inteiras) | PR \#3 |
| 2\. Grades de formato | [`shape_grid()`](https://provetelab.org/treeshapesim/reference/shape_grid.md): cruza níveis de γ e Colless e dá β e `c` de cada célula | PR \#3 |
| 3\. Várias árvores por chamada | `nsim` em [`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md), que devolve `multiPhylo` | PR \#3 |
| 4\. Profundidade da raiz | `root_age` em [`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md) e [`sim_tree_target()`](https://provetelab.org/treeshapesim/reference/sim_tree_target.md) | PR \#3 |
| 5\. Desempenho | **Encerrado.** Benchmark: código compilado não compensa (ver abaixo). Antes: [`assign_ages()`](https://provetelab.org/treeshapesim/reference/assign_ages.md) ficou sem [`ape::prop.part()`](https://rdrr.io/pkg/ape/man/boot.phylo.html) (ganho pequeno, árvores idênticas). O gargalo é [`rtopology_beta()`](https://provetelab.org/treeshapesim/reference/rtopology_beta.md), O(n²) em árvores desbalanceadas; pré-calcular os `lgamma` mudava algumas árvores por arredondamento e foi revertido | PR \#3 |
| 6\. Vignette | [`vignette("treeshapesim")`](https://provetelab.org/treeshapesim/articles/treeshapesim.md): separação dos processos, calibração, alvos e grade | PR \#3 |
| 7\. Citação | `inst/CITATION` e README com o DOI de conceito do Zenodo | PR \#3 e commit seguinte |
| 11\. [`tree_shape()`](https://provetelab.org/treeshapesim/reference/tree_shape.md) | γ e Colless normalizado (Yule) em R puro, iguais a [`ape::gammaStat()`](https://rdrr.io/pkg/ape/man/gammaStat.html) e `treestats::colless(normalization = "yule")`; o treestats continua em Suggests | PR \#3 |
| 13\. Validar [`assign_ages()`](https://provetelab.org/treeshapesim/reference/assign_ages.md) | intervalos positivos e finitos; árvore binária e enraizada | PR \#3 |
| Comparação com outros pacotes | tabela no README (página inicial) com TreeSim e apTreeshape | commit seguinte |
| Release 0.2.0 | versão 0.2.0 no DESCRIPTION e no NEWS; release no GitHub; DOI da versão 10.5281/zenodo.23666434 (conceito 10.5281/zenodo.23647260) | release 0.2.0 |
| Release 0.2.1 | teste do [`fit_shape()`](https://provetelab.org/treeshapesim/reference/fit_shape.md) com árvore fixa (as topologias simuladas diferem entre arm64 e x86-64 para a mesma semente) e nota de reprodutibilidade em [`?sim_tree_shape`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md); CI verde em macOS, Windows e Ubuntu (devel, release e oldrel); DOI da versão 10.5281/zenodo.23667814 | release 0.2.1 |
| [`fit_shape()`](https://provetelab.org/treeshapesim/reference/fit_shape.md) | β (β-splitting, a partir da topologia) e `c` (modelo DDX do laser, com c = −x) de árvores empíricas, com bootstrap paramétrico, detecção de limites e tolerância de ultrametricidade; conferido contra o apTreeshape e o laser arquivados (diferenças ≤ 3e-5) | PR \#4 |
| Validação dos estimadores | `data-raw/validate_estimators.R`: o gerador segue o β-splitting; o DDX recupera `c` com viés ≤ 0,06 | commit `74fb2a2` |
| 9\. Taxa de especiação definida pelo usuário | argumento `rate` em [`intervals_lineage_rate()`](https://provetelab.org/treeshapesim/reference/intervals_lineage_rate.md) e [`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md): função de k com a taxa por linhagem, validada em 2:n; `rate = NULL` mantém as árvores; `rate` com `c` diferente de 0 dá erro; [`calibrate_shape()`](https://provetelab.org/treeshapesim/reference/calibrate_shape.md), [`shape_params()`](https://provetelab.org/treeshapesim/reference/shape_params.md) e [`sim_tree_target()`](https://provetelab.org/treeshapesim/reference/sim_tree_target.md) continuam específicos do `c` | versão de desenvolvimento |
| 8\. Fração amostrada | `sampling_fraction` em [`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md): árvore completa de `round(n / sampling_fraction)` pontas, das quais `n` ficam ao acaso, com profundidade `root_age`; [`sim_tree_target()`](https://provetelab.org/treeshapesim/reference/sim_tree_target.md) e a calibração não tratam disso | versão de desenvolvimento |
| 10\. Parâmetros na árvore | atributo `params` (`beta`, `c`, `root_age`, `sampling_fraction`) e atributo `rate`, como em [`sim_tree_target()`](https://provetelab.org/treeshapesim/reference/sim_tree_target.md) | versão de desenvolvimento |
| 12\. Aviso para `c` extremo | aviso quando há ramos de comprimento zero (c = −20 com 100 pontas já gera); erro claro quando as taxas estouram ou zeram (c = ±400) | versão de desenvolvimento |
| Integração com o comdynamite | seção no README: árvores do treeshapesim como species pools do comdynamite | versão de desenvolvimento |
| 17\. Cobertura de testes | workflow `test-coverage.yaml` (`covr` e Codecov), `codecov.yml` e badge no README; repositório ativo no Codecov, com o segredo `CODECOV_TOKEN` | versão de desenvolvimento |
| [`fit_shape()`](https://provetelab.org/treeshapesim/reference/fit_shape.md) com amostragem incompleta | `sampling_fraction` e `correction_reps`: `c` corrigido por simulação (curva monótona invertida), também no bootstrap; β não é corrigido | versão de desenvolvimento |
| r-universe | pacote listado em `diogoprov/diogoprov.r-universe.dev` (`packages.json`, branch `main`); a versão 0.2.1 está publicada em diogoprov.r-universe.dev, com badge e instrução de instalação no README | `packages.json` |
| 14\. Documentação precisa | a independência de γ e Colless em relação a `c` e β agora é descrita como exata (para um dado n) | PR \#3 |

## A fazer

### Alta prioridade

### Média prioridade

- **Calibração em cache.** Guardar curvas pré-calculadas para tamanhos
  comuns (por exemplo, 50, 100 e 200 pontas), para evitar a calibração a
  cada sessão. Exige uma rodada de simulação só para isso.
- **Integração com o apTreeshape na vignette.**
  [`assign_ages()`](https://provetelab.org/treeshapesim/reference/assign_ages.md)
  aceita qualquer topologia binária; um parágrafo mostrando topologias
  do apTreeshape (convertidas com `as.phylo()`) com os tempos de
  ramificação do treeshapesim.

### Baixa prioridade

- **Conferir o r-universe depois do merge.** Não há submissão a fazer: o
  pacote já está lá e acompanha a `main`. Falta confirmar que a versão
  de desenvolvimento (hoje 0.2.1.9000) compilou sem erro no painel de
  builds.

- **Calibração com amostragem incompleta.** Só se houver uso para grades
  de γ com amostragem: recalcular a curva γ(c) em árvores podadas.
  [`sim_tree_target()`](https://provetelab.org/treeshapesim/reference/sim_tree_target.md)
  não deve ser estendido, porque a poda quebra a independência em que a
  rejeição separada se apoia (o γ passa a depender um pouco da
  topologia).

- **Manter a tabela comparativa do README** atualizada quando TreeSim ou
  apTreeshape mudarem (hoje cita TreeSim 2.4 e apTreeshape 1.5-0.1, este
  arquivado no CRAN em 2023).

- **16. CRAN.** Faz mais sentido depois do aceite do artigo; exige
  `cran-comments.md`, as referências no campo Description no formato
  `<doi:...>` e checks no win-builder e no R-hub.

- **18. Revisão no rOpenSci.**

Retirado: **15. Outros modelos de topologia.** O apTreeshape já oferece
Yule, PDA, Aldous’ branching, biased, funções de partição do usuário e o
modelo de Maliet et al.; em vez de duplicar, basta a integração descrita
acima.

## Desempenho: código compilado não compensa (benchmark de 2026-10-26)

`data-raw/benchmark_topology.R`, máquina ociosa (carga 2,4), Apple
Silicon, R 4.6.1, mediana de 3 repetições; resultados em
`data-raw/benchmark/`. O protótipo em C++
(`data-raw/prototype/rtopology_beta.cpp`) reproduz exatamente as árvores
e o estado do gerador da versão em R (usa `unif_rand`, `lgammafn` e
`revsort` do R).

| n = 20.000 | topologia em R | em C++ | ganho | [`assign_ages()`](https://provetelab.org/treeshapesim/reference/assign_ages.md) |
|----|---:|---:|---:|---:|
| β = −1,9 | 5,4 s | 4,5 s | 1,2× | 0,4 s |
| β = −1 | 0,14 s | 0,05 s | 2,8× | 2,1 s |
| β = 0 | 0,11 s | 0,02 s | 4,6× | 2,4 s |
| β = 2 | 0,11 s | 0,02 s | 5,3× | 2,4 s |

- Onde a topologia é cara (árvores muito desbalanceadas), o custo é o
  cálculo de O(n²) termos `lgamma`, e o C++ só ganha 1,2 a 2×.
  Reaproveitar esses termos mudaria algumas árvores por arredondamento,
  o que já foi descartado.
- Onde o C++ ganha 5×, a topologia já leva cerca de 0,1 s; o tempo total
  passa a ser dominado por
  [`assign_ages()`](https://provetelab.org/treeshapesim/reference/assign_ages.md)
  (cerca de 2,4 s com 20.000 pontas em árvores equilibradas), que também
  precisa repetir o [`sample()`](https://rdrr.io/r/base/sample.html) do
  R para manter as árvores.
- Conclusão: o pacote fica em R puro. Mesmo o pior caso testado (20.000
  pontas, β = −1,9) leva cerca de 6 s por árvore. O protótipo fica em
  `data-raw/` para o caso de alguém precisar de árvores muito maiores;
  nesse caso, C ou C++ via API do R seria o caminho (Julia e Rust não
  conseguem reproduzir as árvores sem reimplementar o gerador e o
  [`sample()`](https://rdrr.io/r/base/sample.html) do R).

## Validação dos estimadores (rodada em 2026-10-26)

`data-raw/validate_estimators.R`, 200 árvores por condição, n = 50, 100
e 200; resultados em `data-raw/validation/` (`summary.csv`). Rodou em 4
min, num núcleo, depois que a máquina liberou. O apTreeshape 1.5-0.1 foi
instalado do arquivo do CRAN numa biblioteca separada.

- **β por `apTreeshape::maxlik.betasplit()`.** Para β ≤ −1, praticamente
  sem viés (≤ 0,15 com n = 50; ≈ 0 com n = 200). Para β ≥ 0, viés
  positivo e grande variância em árvores pequenas (n = 50, β = 0: média
  0,44, REQM 1,38; com β = 3, estimativas no limite superior da busca,
  10), que diminuem com n (β = 0: viés 0,06 com n = 200). É o
  comportamento de um estimador consistente aplicado ao modelo certo:
  **o gerador de topologias segue o β-splitting.** A imprecisão para
  árvores equilibradas é do estimador.
- **`c` pelo DDX.** Viés pequeno e negativo (−0,04 a −0,06 com n = 50;
  −0,01 a −0,02 com n = 200) e REQM de ≈ 0,2 (n = 50) a ≈ 0,08 (n =
  200). O viés vem de o pacote sortear o último intervalo da mesma lei,
  que o DDX não trata como evento; pode ser corrigido no
  [`fit_shape()`](https://provetelab.org/treeshapesim/reference/fit_shape.md)
  ou removido pelo bootstrap paramétrico.

## Correção de `c` para amostragem incompleta (rodada em 2026-11-05)

`data-raw/validate_sampling.R`: n = 100, β = 0, 40 árvores por condição,
`correction_reps = 50`; resultados em
`data-raw/validation/sampling_summary.csv`. Levou 12 min em 4 núcleos,
cerca de 4 s por ajuste.

| fração | `c` verdadeiro | viés sem correção | viés corrigido | REQM sem correção | REQM corrigida |
|----|---:|---:|---:|---:|---:|
| 0,75 | 0 | −0,08 | 0,03 | 0,14 | 0,10 |
| 0,75 | −1 | −0,16 | 0,03 | 0,19 | 0,10 |
| 0,5 | 0 | −0,22 | 0,00 | 0,25 | 0,09 |
| 0,5 | −1 | −0,42 | 0,00 | 0,44 | 0,09 |
| 0,25 | 0 | −0,40 | 0,00 | 0,41 | 0,06 |
| 0,25 | −1 | −0,78 | −0,01 | 0,78 | 0,07 |

Não foram testados outros tamanhos, β ≠ 0, nem a cobertura dos
intervalos do bootstrap. A correção supõe pontas ausentes ao acaso e
fração conhecida.

## O que a literatura diz sobre estimar β e `c`

Busca feita em 2026-10-26. Referências conferidas nas fontes indicadas.

- **O estimador de `c` já existe: é o modelo DDX.** No pacote laser
  (Rabosky 2006a, <doi:10.1177/117693430600200024>), a taxa de
  especiação por linhagem no modelo DDX é r·k^(−x) com k linhagens
  (conferido no código-fonte de `laser::DDX()`). No treeshapesim, o
  intervalo com k linhagens tem taxa total k^(1 + c). Portanto, **c =
  −x**, e estimar `c` por máxima verossimilhança é ajustar o DDX. O
  laser saiu do CRAN, mas o modelo continua sendo usado. Um
  [`fit_shape()`](https://provetelab.org/treeshapesim/reference/fit_shape.md)
  deve citar o DDX e o laser, e não apresentar a estimação como
  novidade. Em dados com o último intervalo sorteado da mesma lei (como
  em
  [`intervals_lineage_rate()`](https://provetelab.org/treeshapesim/reference/intervals_lineage_rate.md)),
  um pequeno viés é possível; a validação agendada mede isso.
- **Cuidados de interpretação.**
  - Escolher entre modelos de diversificação pelo menor AIC infla muito
    o erro tipo I (Rabosky 2006b, <doi:10.1554/05-424.1>).
  - Em modelos de diversidade-dependência, as estimativas podem ser
    enviesadas, e o teste da razão de verossimilhança e o AIC não são
    confiáveis; a recomendação é bootstrap paramétrico (Etienne et
    al. 2016, <doi:10.1111/2041-210x.12565>). O simulador do pacote
    serve exatamente para esse bootstrap.
  - Padrões de ramificação não distinguem diversidade-dependência de
    dependência do tempo (Pannetier et al. 2020,
    <doi:10.1111/evo.14124>). Logo, um `c` estimado descreve a forma da
    curva de acúmulo de linhagens, não o mecanismo. Isso precisa estar
    na documentação.
- **O estimador de β também existe.** O β do modelo de β-splitting é
  estimado por máxima verossimilhança a partir da topologia e foi
  aplicado a um grande conjunto de filogenias publicadas por Blum &
  François (2006, <doi:10.1080/10635150600889625>); está implementado em
  `apTreeshape::maxlik.betasplit()`. Como o apTreeshape foi arquivado no
  CRAN em 2023, o treeshapesim não pode depender dele; a
  verossimilhança, que é simples, deve ser reimplementada e citada.
- **Validade da validação cruzada.** Como `maxlik.betasplit()`
  implementa a verossimilhança do mesmo modelo de Aldous (1996),
  recuperar o β verdadeiro em árvores de
  [`rtopology_beta()`](https://provetelab.org/treeshapesim/reference/rtopology_beta.md)
  é um teste legítimo de que o gerador está correto. No caso do `c`,
  comparar com o DDX é um teste de consistência do gerador de intervalos
  e, ao mesmo tempo, mede o viés do estimador nas árvores do pacote.
- **Conclusão.** As duas estimações são válidas e já publicadas. O valor
  que o treeshapesim acrescenta é juntá-las numa função (β e `c` de uma
  árvore empírica, na mesma escala da simulação) e oferecer o bootstrap
  paramétrico com o próprio gerador, não um método novo.

## Ligação com o artigo

- O `pcps-sim` usa `R/trees.R`, que gera as mesmas árvores que
  [`sim_tree_shape()`](https://provetelab.org/treeshapesim/reference/sim_tree_shape.md)
  com a mesma semente; o teste com valores congelados protege isso.
- Depois do aceite: citar o treeshapesim nos Métodos ou no Apêndice S1,
  e só então ligar o pacote ao artigo no README e no site (anonimato da
  revisão).

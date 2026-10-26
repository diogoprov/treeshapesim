# treeshapesim: estado e próximos passos

Atualizado em 2026-10-26, depois do merge do PR #3, da comparação com TreeSim e
apTreeshape (README) e de uma busca na literatura (Consensus, Undermind e
Scholar Gateway) sobre estimar os parâmetros de árvores empíricas. Versão de
desenvolvimento 0.1.0.9000.

## Feito

| Item | O que entrou | Onde |
|---|---|---|
| Validação e robustez | `n`, `beta` e `c` validados (`c = Inf` dava ramos `NaN`); topologia construída sem recursão (n = 3.000 com β = −1,99 deixou de falhar), com as mesmas árvores para a mesma semente | PR #2 |
| R CMD check na CI | workflow `R-CMD-check.yaml` (r-lib/actions v2) e badge no README | PR #2 |
| 1. Simular para valores-alvo | `calibrate_shape()` (curvas γ(c) e Colless(β), monótonas por regressão isotônica), `shape_params()` (inversão) e `sim_tree_target()` (rejeição separada de topologias e intervalos, equivalente à rejeição de árvores inteiras) | PR #3 |
| 2. Grades de formato | `shape_grid()`: cruza níveis de γ e Colless e dá β e `c` de cada célula | PR #3 |
| 3. Várias árvores por chamada | `nsim` em `sim_tree_shape()`, que devolve `multiPhylo` | PR #3 |
| 4. Profundidade da raiz | `root_age` em `sim_tree_shape()` e `sim_tree_target()` | PR #3 |
| 5. Desempenho | **Parcial.** `assign_ages()` ficou sem `ape::prop.part()` (ganho pequeno, árvores idênticas). O gargalo é `rtopology_beta()`, O(n²) em árvores desbalanceadas; pré-calcular os `lgamma` mudava algumas árvores por arredondamento e foi revertido | PR #3 |
| 6. Vignette | `vignette("treeshapesim")`: separação dos processos, calibração, alvos e grade | PR #3 |
| 7. Citação | `inst/CITATION` e README com o DOI de conceito do Zenodo | PR #3 e commit seguinte |
| 11. `tree_shape()` | γ e Colless normalizado (Yule) em R puro, iguais a `ape::gammaStat()` e `treestats::colless(normalization = "yule")`; o treestats continua em Suggests | PR #3 |
| 13. Validar `assign_ages()` | intervalos positivos e finitos; árvore binária e enraizada | PR #3 |
| Comparação com outros pacotes | tabela no README (página inicial) com TreeSim e apTreeshape | commit seguinte |
| 14. Documentação precisa | a independência de γ e Colless em relação a `c` e β agora é descrita como exata (para um dado n) | PR #3 |

## A fazer

### Alta prioridade

- **Release 0.2.0.** O DOI de conceito do Zenodo (10.5281/zenodo.23647260) já
  está no `inst/CITATION` e no README; a release 0.1.0 (DOI de versão
  10.5281/zenodo.23647261) é anterior aos PRs #2 e #3. Trocar 0.1.0.9000 por
  0.2.0 no DESCRIPTION e no NEWS e publicar a release, que gera o DOI da 0.2.0.
- **Desempenho da topologia (resto do item 5).** Só vale se alguém precisar de
  árvores muito grandes e desbalanceadas. Para não mudar as árvores existentes,
  o caminho seria um gerador novo e opcional (por exemplo, em C++ ou com outro
  algoritmo de amostragem), mantendo o atual como padrão.

### Média prioridade

- **9. Taxa de especiação definida pelo usuário.** Plano detalhado abaixo.
- **Estimar β e `c` de árvores empíricas (`fit_shape()`).** Implementar no
  próprio pacote, sem depender do apTreeshape (arquivado no CRAN em 2023): a
  verossimilhança do β-splitting a partir da topologia (Aldous 1996; usada por
  Blum & François 2006) e a do DDX para `c` (Rabosky 2006). Bootstrap
  paramétrico com o próprio simulador para intervalos de confiança, e aviso
  para árvores pequenas e equilibradas, onde o estimador de β é impreciso (ver a
  validação). Testar contra o apTreeshape arquivado e contra o DDX do laser
  arquivado, só no `data-raw/`, fora do pacote.
- **8. Fração amostrada.** Um argumento para podar pontas ao acaso, porque o γ é
  sensível à amostragem incompleta. A extinção (nascimento e morte) sai do
  escopo, porque é o terreno do TreeSim.
- **10. Guardar os parâmetros na árvore.** Gravar `beta`, `c` e `root_age` como
  atributos do `phylo` em `sim_tree_shape()`; `sim_tree_target()` já guarda os
  seus.
- **12. Aviso para `c` extremo.** Valores finitos muito grandes (por exemplo,
  400) geram ramos de comprimento zero, sem aviso.
- **Calibração em cache.** Guardar curvas pré-calculadas para tamanhos comuns
  (por exemplo, 50, 100 e 200 pontas), para evitar a calibração a cada sessão.
  Exige uma rodada de simulação só para isso.
- **Integração com o apTreeshape na vignette.** `assign_ages()` aceita qualquer
  topologia binária; um parágrafo mostrando topologias do apTreeshape
  (convertidas com `as.phylo()`) com os tempos de ramificação do treeshapesim.

### Baixa prioridade

- **Manter a tabela comparativa do README** atualizada quando TreeSim ou
  apTreeshape mudarem (hoje cita TreeSim 2.4 e apTreeshape 1.5-0.1, este
  arquivado no CRAN em 2023).
- **16. CRAN.** Faz mais sentido depois do aceite do artigo; exige
  `cran-comments.md`, as referências no campo Description no formato
  `<doi:...>` e checks no win-builder e no R-hub.
- **17. Cobertura de testes** (`covr` e badge).
- **18. Revisão no rOpenSci.**
- **19. Partes em C++**, se o desempenho continuar sendo um problema.

Retirado: **15. Outros modelos de topologia.** O apTreeshape já oferece Yule,
PDA, Aldous' branching, biased, funções de partição do usuário e o modelo de
Maliet et al.; em vez de duplicar, basta a integração descrita acima.

## Validação dos estimadores (rodada em 2026-10-26)

`data-raw/validate_estimators.R`, 200 árvores por condição, n = 50, 100 e 200;
resultados em `data-raw/validation/` (`summary.csv`). Rodou em 4 min, num núcleo,
depois que a máquina liberou. O apTreeshape 1.5-0.1 foi instalado do arquivo do
CRAN numa biblioteca separada.

- **β por `apTreeshape::maxlik.betasplit()`.** Para β ≤ −1, praticamente sem viés
  (≤ 0,15 com n = 50; ≈ 0 com n = 200). Para β ≥ 0, viés positivo e grande
  variância em árvores pequenas (n = 50, β = 0: média 0,44, REQM 1,38; com β = 3,
  estimativas no limite superior da busca, 10), que diminuem com n (β = 0:
  viés 0,06 com n = 200). É o comportamento de um estimador consistente aplicado
  ao modelo certo: **o gerador de topologias segue o β-splitting.** A
  imprecisão para árvores equilibradas é do estimador.
- **`c` pelo DDX.** Viés pequeno e negativo (−0,04 a −0,06 com n = 50; −0,01 a
  −0,02 com n = 200) e REQM de ≈ 0,2 (n = 50) a ≈ 0,08 (n = 200). O viés vem de
  o pacote sortear o último intervalo da mesma lei, que o DDX não trata como
  evento; pode ser corrigido no `fit_shape()` ou removido pelo bootstrap
  paramétrico.

## O que a literatura diz sobre estimar β e `c`

Busca feita em 2026-10-26. Referências conferidas nas fontes indicadas.

- **O estimador de `c` já existe: é o modelo DDX.** No pacote laser (Rabosky
  2006a, doi:10.1177/117693430600200024), a taxa de especiação por linhagem no
  modelo DDX é r·k^(−x) com k linhagens (conferido no código-fonte de
  `laser::DDX()`). No treeshapesim, o intervalo com k linhagens tem taxa total
  k^(1 + c). Portanto, **c = −x**, e estimar `c` por máxima verossimilhança é
  ajustar o DDX. O laser saiu do CRAN, mas o modelo continua sendo usado. Um
  `fit_shape()` deve citar o DDX e o laser, e não apresentar a estimação como
  novidade. Em dados com o último intervalo sorteado da mesma lei (como em
  `intervals_lineage_rate()`), um pequeno viés é possível; a validação agendada
  mede isso.
- **Cuidados de interpretação.**
  - Escolher entre modelos de diversificação pelo menor AIC infla muito o erro
    tipo I (Rabosky 2006b, doi:10.1554/05-424.1).
  - Em modelos de diversidade-dependência, as estimativas podem ser enviesadas,
    e o teste da razão de verossimilhança e o AIC não são confiáveis; a
    recomendação é bootstrap paramétrico (Etienne et al. 2016,
    doi:10.1111/2041-210x.12565). O simulador do pacote serve exatamente para
    esse bootstrap.
  - Padrões de ramificação não distinguem diversidade-dependência de
    dependência do tempo (Pannetier et al. 2020, doi:10.1111/evo.14124). Logo,
    um `c` estimado descreve a forma da curva de acúmulo de linhagens, não o
    mecanismo. Isso precisa estar na documentação.
- **O estimador de β também existe.** O β do modelo de β-splitting é estimado por
  máxima verossimilhança a partir da topologia e foi aplicado a um grande
  conjunto de filogenias publicadas por Blum & François (2006,
  doi:10.1080/10635150600889625); está implementado em
  `apTreeshape::maxlik.betasplit()`. Como o apTreeshape foi arquivado no CRAN
  em 2023, o treeshapesim não pode depender dele; a verossimilhança, que é
  simples, deve ser reimplementada e citada.
- **Validade da validação cruzada.** Como `maxlik.betasplit()` implementa a
  verossimilhança do mesmo modelo de Aldous (1996), recuperar o β verdadeiro em
  árvores de `rtopology_beta()` é um teste legítimo de que o gerador está
  correto. No caso do `c`, comparar com o DDX é um teste de consistência do
  gerador de intervalos e, ao mesmo tempo, mede o viés do estimador nas árvores
  do pacote.
- **Conclusão.** As duas estimações são válidas e já publicadas. O valor que o
  treeshapesim acrescenta é juntá-las numa função (β e `c` de uma árvore
  empírica, na mesma escala da simulação) e oferecer o bootstrap paramétrico com
  o próprio gerador, não um método novo.

## Plano do item 9: taxa de especiação definida pelo usuário

**Objetivo.** Permitir qualquer dependência da taxa de especiação por linhagem
em relação ao número de linhagens, `rate(k)`, mantendo os intervalos
independentes da topologia (o que separa o pacote do TreeSim).

**Interface.**

- `intervals_lineage_rate(n, c = 0, rate = NULL)` e
  `sim_tree_shape(n, beta = 0, c = 0, nsim = 1, root_age = 1, rate = NULL)`.
- `rate` é uma função vetorizada de k que devolve a taxa de especiação por
  linhagem com k linhagens. O intervalo com k linhagens é então exponencial com
  taxa k·rate(k).
- Com `rate = NULL`, o código atual roda sem mudança, então os padrões continuam
  gerando as mesmas árvores, e o teste com valores congelados continua valendo.
- `rate` e `c` diferente de 0 ao mesmo tempo dão erro (`c` é o caso particular
  rate(k) = k^c).

**Validação de `rate`.** Uma função que, avaliada em 2:n, devolva um vetor do
mesmo comprimento, numérico, finito e positivo; caso contrário, erro com
mensagem clara (por exemplo, um modelo logístico com K ≤ n dá taxa ≤ 0).

**Documentação.**

- Só a forma de rate(k) importa, não a escala: multiplicar a taxa por uma
  constante não muda a árvore, porque a profundidade é reescalada para
  `root_age`.
- Exemplos: rate(k) = k^c (o padrão); rate(k) = 1 − k/K, com K > n (logístico,
  como o DDL do laser); uma taxa com degrau num certo número de linhagens.
- A mesma ressalva de Pannetier et al. (2020): a forma de rate(k) não identifica
  o mecanismo.
- `calibrate_shape()`, `shape_params()` e `sim_tree_target()` continuam
  específicos do `c`; com `rate`, a calibração é responsabilidade de quem usa.
  Isso fica dito explicitamente.

**Testes.**

- Com rate(k) = k^c, a distribuição do γ coincide com a de `c` (mesma semente,
  diferença só de arredondamento; tolerância numérica).
- Com `rate = NULL`, as árvores são idênticas às atuais (teste congelado
  existente).
- Erros para funções que não devolvem n − 1 valores positivos e finitos, e para
  `rate` junto com `c` diferente de 0.
- O modelo logístico com K próximo de n dá γ médio negativo (semente fixa, sem
  limiares frágeis).

**Arquivos.** `R/simulate.R` (argumento e checagem), testes em
`tests/testthat/test-simulate.R`, um parágrafo na vignette, NEWS e README (uma
linha). Sem dependência nova.

**Computação.** Nenhuma rodada pesada; só os testes unitários.

**Esforço estimado.** Cerca de 1 a 2 h, incluindo revisão de código.

## Ligação com o artigo

- O `pcps-sim` usa `R/trees.R`, que gera as mesmas árvores que
  `sim_tree_shape()` com a mesma semente; o teste com valores congelados protege
  isso.
- Depois do aceite: citar o treeshapesim nos Métodos ou no Apêndice S1, e só
  então ligar o pacote ao artigo no README e no site (anonimato da revisão).

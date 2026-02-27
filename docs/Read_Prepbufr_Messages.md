# Mensagens do código `read_prepbufr.f90`

Este documento lista todas as mensagens (`WRITE`, `print*`, etc.) encontradas no arquivo **read_prepbufr.f90**, com suas respectivas linhas aproximadas, significado e contexto de uso.

---

## 1. Erros e interrupções fatais

| Linha              | Mensagem                                          | Significado                                             | Motivo                                                                 |
| ------------------ | ------------------------------------------------- | ------------------------------------------------------- | ---------------------------------------------------------------------- |
| **518**            | `illegal obs type in READ_PREPBUFR`               | Tipo de observação não reconhecido.                     | O código caiu no `else` final da seleção de tipos e chama `stop2(94)`. |
| **637**            | `READ_PREPBUFR: messages exceed maximum`          | Número de mensagens BUFR ultrapassou o limite definido. | Proteção de tamanho de vetor; aciona `stop2(50)`.                      |
| **644**            | `READ_PREPBUFR: reports exceed maximum`           | Número de reports/subconjuntos excedeu o limite.        | Evita overflow de arrays; `stop2(50)`.                                 |
| **1584**           | `***ERROR*** tail number exceeds maximum`         | Estouro no número de identificadores de aeronaves.      | Controle de arrays de "tail numbers"; aciona `stop2(341)`.             |
| **2957**           | `PREPBUFR: mix up in read_prepbufr ,ndata,icount` | Inconsistência entre o total esperado e o contado.      | Checagem final de integridade; `stop2(50)`.                            |
| **3287**, **3350** | `error in SONDE_EXT levs > 255`                   | Número de níveis de sondagem ultrapassou 255.           | Proteção contra excesso de níveis no array.                            |

---

## 2. Avisos e inconsistências não fatais

| Linha         | Mensagem                                         | Significado                                            | Motivo                                                                        |
| ------------- | ------------------------------------------------ | ------------------------------------------------------ | ----------------------------------------------------------------------------- |
| **599**       | `no matching obstype found in obsinfo`           | Tipo de observação não encontrado na tabela `obsinfo`. | A busca retornou `ntmatch==0`; o tipo é ignorado.                             |
| **1173–1175** | `WARNING!! psob: cannot find subtyep...`         | Subtipo ausente na tabela de erros para `psob`.        | Usa coluna 0 da tabela de erros.                                              |
| **1225–1227** | `WARNING!! tob: cannot find subtyep...`          | Subtipo ausente para `tob`.                            | Usa coluna 0 da tabela.                                                       |
| **1276–1279** | `WARNING!! qob: cannot find subtyep...`          | Subtipo ausente para `qob`.                            | Usa coluna 0 da tabela.                                                       |
| **1382–1384** | `WARNING!! pwob: cannot find subtyep...`         | Subtipo ausente para `pwob`.                           | Usa coluna 0 da tabela.                                                       |
| **1525–1527** | `***WARNING*** invalid pressure pob=...`         | Pressão POB inválida (zero, negativa ou ausente).      | Corrige o valor para `bmiss` e prossegue.                                     |
| **1965**      | `warning: obs pressure is too small:`            | Valor de pressão muito pequeno.                        | Verificação adicional de sanidade.                                            |
| **2050–2051** | `***WARNING*** ndata > maxobs for ...`           | Número de observações excede `maxobs`.                 | Corrige `ndata = maxobs` para evitar overflow.                                |
| **3014–3036** | `WARING convinfo file settings are not right...` | Parâmetros conflitantes na tabela `convinfo`.          | Usa a primeira entrada válida.                                                |
| **3092**      | `something is wrong,lat,lon,prest=`              | Valor fora de faixa para variável.                     | Checagem de limite antes de preencher buckets.                                |

---

## 3. Mensagens informativas e de depuração

| Linha    | Mensagem                                 | Significado                                       | Motivo                             |
| -------- | ---------------------------------------- | ------------------------------------------------- | ---------------------------------- |
| **686**  | `new vad flag::`                         | Log de ativação de flag VAD wind.                 | Indicador de método alternativo.   |
| **745**  | `blacklist station ...`                  | Estação bloqueada para o tipo indicado.           | Em lista negra (blacklist).        |
| **800**  | `no messages/reports`                    | Nenhum dado válido encontrado.                    | Arquivo vazio ou inconsistente.    |
| **803**  | `messages/reports = ... ntread = ...`    | Contagem final de mensagens e reports.            | Log de resumo.                     |
| **806**  | `time offset is ... hours.`              | Diferença horária do arquivo vs. ciclo.           | Resultado de `time_4dvar`.         |
| **906**  | `at line 779: obstype, ictype, rmesh...` | Dump de parâmetros da `convinfo`.                 | Diagnóstico de configurações.      |
| **2168** | `sliu diffuu,vv::`                       | Mensagem de debug de componentes de vento.        | Usado em testes.                   |
| **3049** | `dentrip,pmesh,rmesh,ndata=`             | Log de parâmetros de *thinning* e total de dados. | Diagnóstico.                       |
| **3166** | `ntest,disterrmax=`                      | Log de controle de QC espacial.                   | Mostra erros máximos de distância. |
| **3168** | `nvtest,vdisterrmax=`                    | Log adicional de QC espacial.                     | Complementa o anterior.            |
| **3171** | `closbf()`                               | Unidade BUFR fechada.                             | Log final.                         |

---

## 4. Resumo

O código `read_prepbufr.f90` possui três categorias principais de mensagens:

* **Fatais:** Encerram o processamento quando há erro estrutural (overflow, tipo desconhecido, inconsistência grave).
* **Avisos:** Corrigem ou ignoram situações anômalas (pressão inválida, subtipo ausente, duplicatas em tabelas).
* **Informativas:** Servem para controle, monitoramento e depuração do fluxo de leitura e filtragem de observações.

Essas mensagens são essenciais para interpretação de logs de execução do GSI/OBSPROC e devem ser mantidas no arquivo de log para auditoria e diagnóstico.

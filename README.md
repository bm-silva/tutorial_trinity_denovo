# Tutorial de montagem de transcritos *de novo* utilizando o Trinity

## Introdução

Nesta aula utilizaremos dados de RNA-Seq correspondentes ao fungo *Schizosaccharomyces pombe*, contendo reads de 76 pares de base *paired-end* correspondentes a duas amostras: Sp_log (crescimento logarítmico) e Sp_plat (fase de platô). Dentro da pasta contendo essas reads, é possível encontrar os arquivos com final **left.fq** e **right.fq** no formato FASTQ formatado no estilo do sequenciador Illumina. Estes dados e tutorial foram baseados em um *workshop* dos criadores do Trinity (https://github.com/trinityrnaseq/RNASeq_Trinity_Tuxedo_Workshop/wiki/Trinity-De-novo-Transcriptome-Assembly-Workshop)
Para este tutorial, é necessário ter instalado os programas: Trinity, Salmon e edgeR. Além de utilizar um sistema Linux.
Caso precise, eu disponibilizei no meu repositório Docker, uma imagem contendo todos os programas necessários: https://hub.docker.com/r/brunomsilva/transcriptainer
A pasta de instalação do Trinity é aqui referida como **$TRINITY_HOME**.

Na pasta principal deste tutorial, estão disponíveis os seguintes arquivos:

* **rnaseq_data**: pasta com todas as reads que serão utilizadas;
* **uniprot_sprot.pep**: arquivo contendo os fastas de proteínas do banco de dados SwissProt - Uniprot;
* **uniprot_sprot.dmnd**: arquivo de banco de dados criado pelo Diamond a partir dos fastas do Uniprot;
* **Diamond**: o alinhador que iremos utilizar para anotar os transcritos montados;
* **samples.txt**: arquivo que descreve as amostras e as condições utilizadas pelo experimento;
* **filtrar.sh**: um scriptzinho para poder pegar as sequencias dos genes diferencialmente expressos lá no final.

## Montagem dos transcritos

Para poder realizar a montagem dos transcritos, vamos utilizar todas as reads para montar um único arquivo fasta contendo todos os transcritos montados. Isto é feito para futuramente utilizarmos este arquivo para realizar a análise de genes diferenciais. Iremos rodar o Trinity com o seguinte comando:

```
Trinity --seqType fq --SS_lib_type RF --left rnaseq_data/Sp_log.left.fq.gz,rnaseq_data/Sp_plat.left.fq.gz --right rnaseq_data/Sp_log.right.fq.gz,rnaseq_data/Sp_plat.right.fq.gz --CPU 4 --max_memory 2G --output trinity_saida/
```
* **--left** (rnaseq_data/Sp_log.left.fq.gz,rnaseq_data/Sp_plat.left.fq.gz): indicamos o caminho de todas as reads no sentido forward;
* **--right** (rnaseq_data/Sp_log.right.fq.gz,rnaseq_data/Sp_plat.right.fq.gz): indicamos o caminho de todas as reads no sentido reverse; Não se esqueça de que estes arquivos devem ser separados por vírgulas, e não espaçamentos.
* **--SS_lib_type** (RF): seleciona o tipo da sua biblioteca, se é paired-end ou single-end;
* **--CPU** (4): quantidade de threads que iremos utilizar;
* **--max_memory** (2G): quantidade máxima de memória que o Trinity irá utilizar;
* **--output** (trinity_saida/): local no qual o Trinity irá salvar todos os dados relacionados à montagem dos transcritos.

Após completar a montagem, será criado uma pasta chamada **trinity_saida**, tal como o seu output foi descrito. Dentro desta pasta é possível encontrar todos os dados relacionados à montagem dos transcritos.

Ao final, o Trinity irá disponibilizar um arquivo chamado Trinity.fasta, contendo todos os transcritos que foram montados:

```
>TRINITY_DN10_c0_g1_i1 len=2250 path=[1:0-209 2:210-2249]
CTGTCTTTATTCCATCTTTTTTAATACTTTTTCTAAGTTTTCTACTTATTTCCCTACATTTTTTTTTTAAATTTCCCCATCTTTGTTTTGTAA
ACATTTTAAGCTTCAGTTTTCTATATACTGCGTTTGATTTACTTTTGACACATTTATTTATCTTAAATATTTCTTGTTATATATATACGCATC
AGACTGGCTAACTTCTCGTTTGTAGGCTTTACAATGCAGAGGCTAATGTGATTGTGCTATATGGAAACTCGATGGACACATCGGCACGATTAT
CAGGTATCGTAGTTCTTTCTGTATCTTCCCCCATTCGGGTAAAGAATATTAAATTACGGCTAAGTGGCCGTTCTTTCGTTTGTTGGGCTGATG
AGTCCCGTCATGCATCGCCTGGTAACAGAATTCGACGTCAAGTTGTACAAATTCTGGATAAAAGCTGGTCTTTCCTTGCTCCCAACGAGTCAG
```

* O identificador **>** é utilizado como todo arquivo fasta, ser um cabeçalho para a sequência;
* **TRINITY_DN94**, como no exemplo, é o nome do transcrito montado;
* **_c0**: indica que este transcrito pertence ao *cluster* 0;
* **_g1**: indica que é o gene 1;
* **_i1**: indica que é a isoforma 1.

Vamos analisar agora os dados relacionados a esta montagem utilizando um *script* do Trinity chamado TrinityStats.pl:

```
TrinityStats.pl trinity_saida/Trinity.fasta
```

```
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  396
Total trinity transcripts:      400
Percent GC: 39.16
########################################
Stats based on ALL transcript contigs:
########################################
        Contig N10: 2589
        Contig N20: 2148
        Contig N30: 1818
        Contig N40: 1594
        Contig N50: 1329
        Median contig length: 548.5
        Average contig: 839.30
        Total assembled bases: 335720
#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################
        Contig N10: 2589
        Contig N20: 2118
        Contig N30: 1818
        Contig N40: 1584
        Contig N50: 1322
        Median contig length: 545.5
        Average contig: 830.32
        Total assembled bases: 328805
``` 

Neste reporte é possível ter acesso à informações como: Número de genes encontrados, número de transcritos, a porcentagem de GC global e a média do tamanho dos contigs.
* Nx: Todos os contigs são ordenados, a partir do seu tamanho, do maior para o menor; o menor contig montado dentre os maiores compreendendo x% dos nucleotídeos totais, representa o seu Nx.

## Abundância e contagem de cada transcrito

Vamos agora quantificar e contar a abundância de cada transcrito e gene que foram montados utilizando o Salmon para cada condição com os seguintes comandos:

```
$TRINITY_HOME/util/align_and_estimate_abundance.pl --seqType fq --left rnaseq_data/Sp_log.left.fq.gz --right rnaseq_data/Sp_log.right.fq.gz --transcripts trinity_saida/Trinity.fasta --est_method salmon --trinity_mode --prep_reference --output_dir salmon_saida/sp_log/
```
```
$TRINITY_HOME/util/align_and_estimate_abundance.pl --seqType fq --left rnaseq_data/Sp_plat.left.fq.gz --right rnaseq_data/Sp_plat.right.fq.gz --transcripts trinity_saida/Trinity.fasta --est_method salmon --trinity_mode --prep_reference --output_dir salmon_saida/sp_plat/
```

* **align_and_estimate_abundance.pl**: script do pacote Trinity que ira realizar todo o pipeline;
* **--seqType** (fq): indicamos qual o formato dos arquivos das reads;
* **--left** (rnaseq_data/Sp_log.left.fq.gz): indicamos o caminho de todas as reads no sentido forward;
* **--right** (rnaseq_data/Sp_log.right.fq.gz): indicamos o caminho de todas as reads no sentido reverse;
* **--transcripts** (trinity_saida/Trinity.fasta): o caminho do arquivo contendo todos os transcritos montados pelo Trinity;
* **--est_method** (salmon): indicamos qual o método/software de estimação iremos utilizar;
* **--trinity_mode**: irá gerar um arquivo que irá mapear todas as isoformas em seus respectivos genes;
* **--prep_reference**: irá criar índices no arquivo dos transcritos para agilizar os processos;
* **--output_dir** (salmon_saida/sp_log/): indicamos qual o local dos arquivos de saída.

o final do processo, irá ser criado uma pasta chamada salmon_saida com duas pastas: Sp_log e Sp_plat, contendo todos os resultados obtidos pelo Salmon.
Dentro de cada pasta de saída do último comando podemos encontrar 2 arquivos importantes, o **quant.sf** e o **quant.sf.genes**:

* **quant.sf**: quantificação da abundância dos transcritos;
* **quant.sf.genes**: quantificação de abundância a nível de genes.

Estes arquivos são separados por tabulações e possuem colunas das quais podemos identificar o nome do transcrito/gene, o tamanho, e o valores de TPM e de suas reads:
```
Name    Length  EffectiveLength TPM     NumReads
TRINITY_DN0_c0_g1_i1    604     335.717 322.040963      6.000
TRINITY_DN0_c0_g2_i1    1838    1569.717        1239.753502     108.000
TRINITY_DN1_c0_g1_i1    563     294.717 54596.288349    892.966
```

Agora vamos comparar os genes e transcritos entre as amostras utilizando um script do pacote Trinity chamado: **abundance_estimates_to_matrix.pl**. 

```
$TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method salmon --out_prefix salmon_saida/salmon salmon_saida/sp_log/quant.sf salmon_saida/sp_plat/quant.sf --gene_trans_map trinity_saida/Trinity.fasta.gene_trans_map --name_sample_by_basedir
```

* **--est_method** (salmon): indicamos qual o método/software de estimação foi utilizado para o cálculo de abundância;
* **--out_prefix** (salmon_saida/salmon): indicamos a parta e o nome dos arquivos de saída;
* **--gene_trans_map** (trinity_saida/Trinity.fasta.gene_trans_map): aquele índice contento cada transcrito pareado com a sua isoforma gerado no passo 5;
* **--name_sample_by_basedir**: utilizado para poder nomear cada condição de acordo com o nome da pasta que as contém.

Na pasta **salmon_saida** podemos ver que apareceram novos arquivos:

* Aqueles que possuem isoform são os resultados a nível dos transcritos e aqueles que possuem gene são a nível de genes;
* Arquivos terminados em .counts.matrix possuem a estimativa da contagem dos fragmentos de RNA-Seq;
* **.TPM.not_cross_norm**: matriz contendo os valores de expressão TPM que não estão normalizadas por amostra;
* **.TMM.EXPR.matrix**: matriz contendo os valores de expressão TMM normalizados.

Nestes arquivos de saída podemos observar que as colunas dos valores foram nomeadas de acordo com a pasta que elas estavam, facilitando o trabalho, nos próximos passos, para identificar qual dado é de qual amostra é condição.

## Análise de Experssão Diferencial

Agora iremos descobrir os genes que estão diferencialmente expressos e montar uns gráficos para melhor visualização destes dados que criamos. Vamos utilizar o programa **edgeR** para este cálculo e outro *script* do Trinity chamado **run_DE_analysis.pl**.

```
$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix salmon_saida/salmon.gene.counts.matrix --method edgeR --output edgeR_saida/ --samples_file samples.txt --dispersion 0.1
```

* **--matrix** (salmon_saida/salmon.gene.counts.matrix): indicamos de qual matriz contém os dados de abundâncias dos transcritos de genes. Neste caso, iremos calcular somente os genes diferenciados e não suas isoformas.
* **--method** (edgeR): indicamos qual programa usar;
* **--samples_file** (samples.txt): um arquivo que descreve as replicatas e condições dos seus experimentos. Já já iremos vê-lo;
* **--dispersion** (0.1): termo utilizado na função binomial negativa do edgeR para estimar a contagem dos genes.

O arquivo **samples.txt** é muito importante nessa etapa uma vez que ele indica qual a condição e suas replicatas para o programa. Ele possui somente duas colunas, com suas linhas separadas por tabulações, e os nomes destas colunas devem estar de acordo com o nomes dados na etapa de estimação de abundância e contagem, dentro do arquivo que possui o final **.counts.matrix**. Exemplo:

```
#condition        #samples
condA             repA1
condA             repA2
condB             repB1
condB             repB2
```

Dentro da pasta **edgeR_saida** já é possível encontrar os arquivos de saída:

* O nome do arquivo segue uma lógica: **{prefix}.sampleA_vs_sampleB.{method}**;
* **.DE.results**: arquivo contendo os resultados das análises, incluindo o fold change e a significância estatística (FDR);
* **.MA_n_Volcano.pdf**: gráficos gerádos pelo *script* do Trinity, um de MA e um Volcano!

Quantos genes estão diferencialmente expressos? Vamos descobrir com o seguinte comando (dentro da pasta **edgeR_saida**):

```
sed '1,1d' salmon.gene.counts.matrix.sp_log_vs_sp_plat.edgeR.DE_results | awk '{ if ($7 <= 0.05) print;}' | wc -l
```

* O comando **sed** é uma ferramenta Unix que nos permite editar arquivos texto. com *'1,1d'* estamos ignorando a primeira linha, o *header*;
* O comando **awk** permite escrever no terminal somente as linhas que contém valores menores que **0.05** na coluna de **FDR** (if ($7 <= 0.05) print);
* Por fim, contamos as linhas com **wc -l**.

## Anotação dos Genes Diferencialmente Expressos

Vamos agora pegar todas as sequências de todos os genes que possuem o FDR <= 0.05 utilizando o *script* **filtrar.sh**. Para isso, digite no terminal:

```
sh filtrar.sh
```

Serão gerados 2 arquivos:

* **DE_genes.txt**: contém o nome de cada gene;
* **DE_genes.fasta**: contém o *header* e a sequência do gene;

Agora vamos utilizar o Diamond para poder identificar a nomenclatura de cada gene para podermos conhecê-los. Para rodar o Diamond, rode o seguinte comando:

```
./Diamond blastx -d uniprot_sprot.dmnd -q DE_genes.fasta -v --outfmt 6 -e 1e-10 --max-target-seqs 1 --outfmt 6 -o blastx.outfmt6 -b4
```

* Estamos rodando o blastx do Diamond;
* **-d** (../data/uniprot_sprot.dmnd): indicamos o caminho do banco de dados que será utilizado para blastar as sequências. Neste caso, estamos utilizando o SwissProt-Uniprot;
* **-q** (DE_genes.fasta): indicamos a query, que neste caso é o arquivo que contém a sequência dos genes que possuem um valor de FDR <= 0.05;
* **-v**: modo verbose, para poder mostrar letrinhas no terminal;
* **-e** (1e-10): valor de e-value;
* **--max-target-seqs** (1): para poder anotar somente o best-hit de cada transcrito;
* **--outfmt** (6): define o modelo do qual será o arquivo de saída, que neste caso, é igual ao formato de saída do programa Blast;
* **-o** (blastx.outfmt6): arquivo de output.

Será gerado um arquivo contendo o resultado do alinhamento, o **blastx.outfmt6**, que possui 12 colunas:

* Coluna 1: Nome da query;
* Coluna 2: Nome da referência utilizada;
* Coluna 3: porcentagem de identidade;
* Coluna 4: tamanho do alinhamento;
* Coluna 5: número de mismatches;
* Coluna 6: número de gaps;
* Coluna 7; posição de início na sequência da query;
* Coluna 8: posição final na sequência da query;
* Coluna 9: posição de início na sequência da referência;
* Coluna 10: posição final na sequência da referência;
* Coluna 11: e-value;
* Coluna 12: bit score.

Para ter acesso a todos os hits, digite no terminal:

```
awk '{print $2}' blastx.outfmt6
```

## Gráficos e Clusterizações

Iremos agora extrair esses genes e gerar os famosos gráficos de *heatmap*. Para isso utilizaremos o **analyze_diff_expr.pl** do Trinity, mas precisamos estar dentro da pasta **edgeR_saida**:
  
```
$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl  -P 1e-3 -C 2 --matrix ../salmon_saida/salmon.gene.TMM.EXPR.matrix --samples ../samples.txt
```

* **-P** (1e-3): define o p-valor a ser utilizado para ser considerado um gene diferencialmente expresso;
* **-C** (2): define o valor do fold change;
* **--matrix** (../salmon_saida/salmon.gene.TMM.EXPR.matrix): indicamos de qual matriz contém os dados de abundâncias dos transcritos de genes gerados no passo 7 contendo os valores de TMM entre as amostras;
* **--samples** (../samples.txt): aquele mesmo arquivo que descreve as replicatas e condições do seus experimentos.

Arquivos importantes criados:

* **.{sampleA}-UP.subset**: features que estão reguladas positiviamentes na amostra A;
* **.{sampleB}-UP.subset**: features que estão reguladas positiviamentes na amostra B;
* **diffExpr.{pvalor}_{valorFC}.matrix.log2.dat**: todas a features encontradas que estão diferencialmente expressas em todas as comparações entre as amostras;
* **diffExpr.{pvalor}_{valorFC}.matrix.log2.sample_cor.dat**: uma matriz de correlação de Pearson para comparações entre as amostras baseada nas features diferencialmente expressas encontradas;
* **diffExpr.{pvalor}_{valorFC}.matrix.log2.centered.genes_vs_samples_heatmap.pdf**: nome muito grande para poder falar que é um gráfico de heatmap clusterizado dos genes diferencialmente expressos em cada replicata de suas amostras;
* **diffExpr.{pvalor}_{valorFC}.matrix.log2.sample_cor_matrix.pdf**: outro heatmap clusterizado mostrando as matrizes de correlações entre as amostras;

Agora vamos terminar as análises de genes diferenciados e fazer uma clusterização dos mesmos para melhor visualização de grupos que com um padrão similar nas amostras.
Dê o seguinte comando, ainda dentro da pasta **edgeR_saida**:

```
$TRINITY_HOME/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl -R diffExpr.P1e-3_C2.matrix.RData --Ktree 3
```

* **-Ktree**: Indica que vamos utilizar valor de k-means para poder fazer esta clusterização. É possível alterar  este valor até achar um que melhor represente o seu conjunto de amostras.










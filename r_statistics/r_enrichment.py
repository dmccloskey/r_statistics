from .r_dependencies import *
from .r_base import r_base
import collections
class r_enrichment(r_base):

    def import_GODB(self,GO_I='GO.db'):
        '''import the GO database of choice into the R workspace
        INPUT:
        GO_I = string, name of the GO library
        '''
        self.library_bioconductor(GO_I);

    def make_topDiffGenes(self,
            topDiffGenes_O = 'topDiffGenes',
            pvalue_I = 0.05,
            ):
        """
        Make the topDiffGenes function to determine differential expressed genes for TopGo
        
        INPUT:
        topDiffGenes_O = string, function to classify genes for topGo
        pvalue_I = float, pvalue threshold
        
        """
        try:
            r_statement = ('%s = function(allScore) {return (allScore < %s)}' % (topDiffGenes_O,pvalue_I));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);
        #try:
        #    r_statement = ('%s = amelia(%s, m=%s)' % (ameliaII_object_O,data_I,n_imputations_I));
        #    ans = robjects.r(r_statement);
        #except Exception as e:
        #    print(e);
        #    exit(-1);

    def make_topGOdata(self,
            topDiffdata_O = 'topGOdata_O',
            ontology = "MF",
            annot = "annFUN.org",
            mapping =  "org.EcK12.eg.db",
            ID = 'symbol',
            geneSel = "topDiffGenes",
            allGenes = "genes_w_pvals",
            gene2GO = "gene_to_go",
            nodeSize = 10,
            ):
        """
        Make the topDiffdata function to determine differential expressed genes for TopGo
        
        INPUT:
        topDiffdata_O = string, topGO data object
        ontology = "MF",
        annot = "annFUN.org",
        geneSel = "topDiffGenes",
        allGenes = "genes_w_pvals",
        nodeSize = 10,
        
        TEST: ans = robjects.r('sum(topDiffGenes(genes_w_pvals))');

        """
        try:
            r_statement = ('%s = new("topGOdata",\
                                ontology = "%s",\
                                annot = %s,\
                                mapping = "%s",\
                                ID = "%s",\
                                allGenes = %s,\
                                geneSel = %s,\
                                nodeSize = %s)'                                
                           %(topDiffdata_O,ontology,
                             annot,
                             mapping,
                             ID,
                             allGenes,
                             geneSel,nodeSize
                             ));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);

    def calculate_enrichment_topGO(self,
            topGOresult_O = 'result',
            topDiffdata_I = 'topGOdata_O',
            algorithm = "classic",
            statistic = "fisher",
            ):
        """
        Run the topGo enrichment test
        
        INPUT:
        topGOresult_O = 'result',
        topDiffdata_I = 'topGOdata_O',
        algorithm = "classic",
        statistic = "fisher",
        
        """
        try:
            r_statement = ('%s = runTest(%s,\
                                algorithm = "%s",\
                                statistic = "%s")'                            
                           %(topGOresult_O,topDiffdata_I,
                             algorithm,
                             statistic
                             ));
            ans = robjects.r(r_statement);
        except Exception as e:
            print(e);
            exit(-1);

    def extract_scores_topGO(self,
            topGOresult_I = 'result',
            topGOscores_O = 'scores',
            ):
        '''Extract out topGO scores
        INPUT:
        topGOresult_I = name of the R workspace variable to contain the output of topGO
        OUTPUT:
        data_O = array of dim 1
        '''
        data_O = None;
        try:
            r_statement = ('%s <- score(%s)' % (topGOscores_O,topGOresult_I));
            ans = robjects.r(r_statement);
            GO_ids = np.array(ans.names);
            pvals = np.array(ans);
            return GO_ids, pvals;
        except Exception as e:
            print(e);
            exit(-1);
        return data_O;

    ##
    def get_go_children(self,go_term, go_term_type):
        """Retrieve all more specific GO children from a starting GO term.
        https://github.com/chapmanb/bcbb/blob/master/stats/diffexp_go_analysis.py
        INPUT:
        go_term = string,
        go_term_type = string,
        OUTPUT:
        children = [], 
        """
        child_map = robjects.r["GO%sCHILDREN" % (go_term_type)]
        children = []
        to_check = [go_term]
        while len(to_check) > 0:
            new_children = []
            for check_term in to_check:
                new_children.extend(list(robjects.r.get(check_term, child_map)))
            new_children = list(set([c for c in new_children if c]))
            children.extend(new_children)
            to_check = new_children
        children = list(set(children))
        return children;

    def parse_go_map_file(self,genes2GOID_I, genes_w_pvals):
        """Parse a GO map file
        https://github.com/chapmanb/bcbb/blob/master/stats/diffexp_go_analysis.py
        INPUT:
        genes2GOID_I = [{}] where go_id = enrichment_class
                                   gene_id = component_name,
        genes_w_pvals = {gene_id:pvalue},
        OUTPUT:
        gene_to_go = {}, 
        go_to_gene = {},
        """
        gene_list = genes_w_pvals.keys()
        gene_to_go = collections.defaultdict(list)
        go_to_gene = collections.defaultdict(list)
        for row in genes2GOID_I:
            gene_id = row['component_name'];
            go_id = row['enrichment_class'];
            if gene_id in gene_list:
                gene_to_go[gene_id].append(go_id)
                go_to_gene[go_id].append(gene_id)
        return dict(gene_to_go), dict(go_to_gene)
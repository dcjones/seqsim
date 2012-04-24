

#include "params.hpp"
#include "yaml-cpp/yaml.h"
#include <fstream>

using namespace std;


/* Fail with dignity when an expected node is missing. */
static const YAML::Node& yaml_get_node(const YAML::Node& node, const char* key)
{
    try {
        return node[key];
    }
    catch (...) {
        fprintf(stderr,
            "seqsim: model specification is missing a %s parameter.\n", key);
        exit(EXIT_FAILURE);
    }
}


/* Try to get a node, returning NULL if it is not present. */
static const YAML::Node* yaml_try_get_node(const YAML::Node& node, const char* key)
{
    try {
        return node.FindValue(key);
    }
    catch (...) {
        return NULL;
    }
}


perturb_params::perturb_params()
    : gene_pr(0.0)
    , gene_mu(0.0)
    , gene_sd(0.0)
    , trans_pr(0.0)
    , trans_mu(0.0)
    , trans_sd(0.0)
{
}


perturb_params::~perturb_params()
{
}


params::params()
    /* some sane defaults */
    : N(1000000)
    , readlen(100)
    , paired(true)
    , strand_specificity(0.5)
    , gene_exp_k(2)
    , size_lower(150)
    , size_upper(170)
    , size_sd(5)
{
    /* TODO: reasonable defaults */

    gene_exp_p.resize(2);
    gene_exp_p[0] = 0.5;
    gene_exp_p[1] = 0.5;

    gene_exp_mu.resize(2);
    gene_exp_mu[0] = 0.1;
    gene_exp_mu[1] = 10.0;

    gene_exp_sd.resize(2);
    gene_exp_sd[0] = 1.0;
    gene_exp_sd[1] = 1.0;
}


params::~params()
{
}

void params::read(const char* fn)
{
    std::ifstream fin(fn);
    YAML::Parser parser(fin);

    unsigned int i;
    YAML::Node node;
    parser.GetNextDocument(node);

    const YAML::Node* value;

    if ((value = yaml_try_get_node(node, "N")))                  *value >> N;
    if ((value = yaml_try_get_node(node, "readlen")))            *value >> readlen;
    if ((value = yaml_try_get_node(node, "paired")))             *value >> paired;
    if ((value = yaml_try_get_node(node, "strand_specificity"))) *value >> strand_specificity;
    if ((value = yaml_try_get_node(node, "size_lower")))         *value >> size_lower;
    if ((value = yaml_try_get_node(node, "size_upper")))         *value >> size_upper;
    if ((value = yaml_try_get_node(node, "size_sd")))            *value >> size_sd;

    if ((value = yaml_try_get_node(node, "gene_exp_k"))) {
        *value >> gene_exp_k;
        gene_exp_p.resize(gene_exp_k);
        gene_exp_mu.resize(gene_exp_k);
        gene_exp_sd.resize(gene_exp_k);
    }

    if ((value = yaml_try_get_node(node, "gene_exp_p"))) {
        for (i = 0; i < gene_exp_k; ++i) {
            (*value)[i] >> gene_exp_p[i];
        }
    }

    if ((value = yaml_try_get_node(node, "gene_exp_mu"))) {
        for (i = 0; i < gene_exp_k; ++i) {
            (*value)[i] >> gene_exp_mu[i];
        }
    }

    if ((value = yaml_try_get_node(node, "gene_exp_sd"))) {
        for (i = 0; i < gene_exp_k; ++i) {
            (*value)[i] >> gene_exp_sd[i];
        }
    }

    if ((value = yaml_try_get_node(node, "trans_exp_alpha"))) {
        *value >> trans_exp_alpha;
    }

    if ((value = yaml_try_get_node(node, "perturbations"))) {
        for (YAML::Iterator j = value->begin(); j != value->end(); ++j) {
            string name;
            j.first() >> name;

            perturb_params P;
            yaml_get_node(j.second(), "gene_pr") >> P.gene_pr;
            yaml_get_node(j.second(), "gene_mu") >> P.gene_mu;
            yaml_get_node(j.second(), "gene_sd") >> P.gene_sd;

            yaml_get_node(j.second(), "trans_pr") >> P.trans_pr;
            yaml_get_node(j.second(), "trans_mu") >> P.trans_mu;
            yaml_get_node(j.second(), "trans_sd") >> P.trans_sd;

            perturbations[name] = P;
        }
    }
}


void params::write(const char* fn)
{
    YAML::Emitter emit;

    emit << YAML::BeginDoc;
    emit << YAML::BeginMap;

    emit << YAML::Key   << "N";
    emit << YAML::Value << (unsigned long long) N;

    emit << YAML::Key   << "paired";
    emit << YAML::Value << paired;

    emit << YAML::Key   << "strand_specificity";
    emit << YAML::Value << strand_specificity;

    emit << YAML::Key   << "gene_exp_k";
    emit << YAML::Value << gene_exp_k;

    emit << YAML::Key    << "gene_exp_p";
    emit << YAML::Value << YAML::BeginSeq;
    for (unsigned int i = 0; i < gene_exp_k; ++i) {
        emit << gene_exp_p[i];
    }
    emit << YAML::Value << YAML::EndSeq;

    emit << YAML::Key    << "gene_exp_mu";
    emit << YAML::Value << YAML::BeginSeq;
    for (unsigned int i = 0; i < gene_exp_k; ++i) {
        emit << gene_exp_mu[i];
    }
    emit << YAML::Value << YAML::EndSeq;

    emit << YAML::Key    << "gene_exp_sd";
    emit << YAML::Value << YAML::BeginSeq;
    for (unsigned int i = 0; i < gene_exp_k; ++i) {
        emit << gene_exp_sd[i];
    }
    emit << YAML::Value << YAML::EndSeq;

    emit << YAML::Key    << "trans_exp_alpha";
    emit << YAML::Value  << trans_exp_alpha;

    emit << YAML::Key   << "size_lower";
    emit << YAML::Value << size_lower;

    emit << YAML::Key   << "size_upper";
    emit << YAML::Value << size_upper;

    emit << YAML::Key   << "size_sd";
    emit << YAML::Value << size_sd;

    emit << YAML::EndMap;
    emit << YAML::EndDoc;

    FILE* f = fopen(fn, "w");
    fwrite(emit.c_str(), sizeof(char), emit.size(), f);
    fclose(f);
}


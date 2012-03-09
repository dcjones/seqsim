

#include "params.hpp"
#include "yaml-cpp/yaml.h"
#include <fstream>


params::params()
    /* some sane defaults */
    : N(1000000)
    , paired(true)
    , strand_specificity(0.5)
    , gene_exp_k(2)
    , size_lower(150)
    , size_upper(170)
    , size_var(5)
{
    /* TODO: reasonable defaults */

    gene_exp_p  = new double [2];
    gene_exp_p[0] = 0.5;
    gene_exp_p[1] = 0.5;

    gene_exp_mu = new double [2];
    gene_exp_mu[0] = 0.1;
    gene_exp_mu[1] = 10.0;

    gene_exp_sd = new double [2];
    gene_exp_sd[0] = 1.0;
    gene_exp_sd[1] = 1.0;

    trans_exp_p  = new double [2];
    trans_exp_p[0] = 0.5;
    trans_exp_p[1] = 0.5;

    trans_exp_mu = new double [2];
    trans_exp_mu[0] = 0.1;
    trans_exp_mu[1] = 10.0;

    trans_exp_sd = new double [2];
    trans_exp_sd[0] = 1.0;
    trans_exp_sd[1] = 1.0;
}


params::~params()
{
    delete [] gene_exp_p;
    delete [] gene_exp_mu;
    delete [] gene_exp_sd;

    delete [] trans_exp_p;
    delete [] trans_exp_mu;
    delete [] trans_exp_sd;
}

void params::read(const char* fn)
{
    std::ifstream fin(fn);
    YAML::Parser parser(fin);

    unsigned int i;
    YAML::Node node;
    parser.GetNextDocument(node);

    const YAML::Node* value;

    if ((value = node.FindValue("N")))                  *value >> N;
    if ((value = node.FindValue("paired")))             *value >> paired;
    if ((value = node.FindValue("strand_specificity"))) *value >> strand_specificity;
    if ((value = node.FindValue("size_lower")))         *value >> size_lower;
    if ((value = node.FindValue("size_upper")))         *value >> size_upper;
    if ((value = node.FindValue("size_var")))           *value >> size_var;

    if ((value = node.FindValue("gene_exp_k"))) {
        *value >> gene_exp_k;
        delete [] gene_exp_mu;
        delete [] gene_exp_sd;
        gene_exp_mu = new double [gene_exp_k];
        gene_exp_sd = new double [gene_exp_k];
    }

    if ((value = node.FindValue("gene_exp_mu"))) {
        for (i = 0; i < gene_exp_k; ++i) {
            (*value)[i] >> gene_exp_mu[i];
        }
    }

    if ((value = node.FindValue("gene_exp_sd"))) {
        for (i = 0; i < gene_exp_k; ++i) {
            (*value)[i] >> gene_exp_sd[i];
        }
    }

    if ((value = node.FindValue("trans_exp_k"))) {
        *value >> trans_exp_k;
        delete [] trans_exp_mu;
        delete [] trans_exp_sd;
        trans_exp_mu = new double [trans_exp_k];
        trans_exp_sd = new double [trans_exp_k];
    }

    if ((value = node.FindValue("trans_exp_mu"))) {
        for (i = 0; i < trans_exp_k; ++i) {
            (*value)[i] >> trans_exp_mu[i];
        }
    }

    if ((value = node.FindValue("trans_exp_sd"))) {
        for (i = 0; i < trans_exp_k; ++i) {
            (*value)[i] >> trans_exp_sd[i];
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

    emit << YAML::Key    << "trans_exp_p";
    emit << YAML::Value << YAML::BeginSeq;
    for (unsigned int i = 0; i < trans_exp_k; ++i) {
        emit << trans_exp_p[i];
    }
    emit << YAML::Value << YAML::EndSeq;

    emit << YAML::Key    << "trans_exp_mu";
    emit << YAML::Value << YAML::BeginSeq;
    for (unsigned int i = 0; i < trans_exp_k; ++i) {
        emit << trans_exp_mu[i];
    }
    emit << YAML::Value << YAML::EndSeq;

    emit << YAML::Key    << "trans_exp_sd";
    emit << YAML::Value << YAML::BeginSeq;
    for (unsigned int i = 0; i < trans_exp_k; ++i) {
        emit << trans_exp_sd[i];
    }
    emit << YAML::Value << YAML::EndSeq;

    emit << YAML::Key   << "size_lower";
    emit << YAML::Value << size_lower;

    emit << YAML::Key   << "size_upper";
    emit << YAML::Value << size_upper;

    emit << YAML::Key   << "size_var";
    emit << YAML::Value << size_var;

    emit << YAML::EndMap;
    emit << YAML::EndDoc;

    FILE* f = fopen(fn, "w");
    fwrite(emit.c_str(), sizeof(char), emit.size(), f);
    fclose(f);
}


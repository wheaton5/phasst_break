#[macro_use]
extern crate clap;
extern crate bio;
extern crate hashbrown;
extern crate phasst_lib;

use bio::io::fasta;
use bio::io::fasta::Record;
use bio::utils::TextSlice;
use std::path::Path;

use phasst_lib::{
    load_assembly_kmers, load_hic, load_hifi, load_linked_read_barcodes, Assembly, Mols, Kmers, KmerMols, Molecules, KmerType,
};
use rayon::prelude::*;


use hashbrown::{HashMap, HashSet};

use clap::App;


use std::fs::File;
use std::io::{BufWriter, Write};
use std::io::BufReader;
use std::io::BufRead;
use std::io::Read;


struct ContigLoci {
    kmers: HashMap<i32, ContigLocus>, // map from kmer id to contig id and position and which allele assembly had
    loci: HashMap<i32, Vec<ContigLocus>>, // map from contig id to number of loci
}

#[derive(Debug, Clone, Copy)]
struct ContigLocus {
    contig_id: i32,
    index: usize,
    position: usize,
    allele: Allele,
    reference: i32,
    alternate: i32,
}

const MIN_ALLELE_FRACTION_HIC: f32 = 0.15;


fn main() {
    println!("Welcome to phasst break!");
    let params = load_params();
    eprintln!("loading kmers");
    let kmers = Kmers::load_kmers(&params.het_kmers);
    //let (_variants, molecules) = load_molecule_kmers(&params.txg_mols, &params.hic_mols, &params.longread_mols, &kmers);
    eprintln!("loading hic kmers");
    let hic_mols = load_hic(Some(&params.hic_mols), &kmers, true); 
    eprintln!("loading assembly kmers");
    let assembly = load_assembly_kmers(&params.assembly_kmers, &params.assembly_fasta, &kmers);
    let allele_fractions = get_allele_fractions(&hic_mols, &kmers); // MAYBE ADD LINKED Molecule AND LONG Molecule to this?
    let bad_alleles = get_bad_alleles(&hic_mols);
    let variant_contig_order: ContigLoci =
        good_assembly_loci(&assembly, &allele_fractions, &bad_alleles);
    let hic_links = gather_hic_links(&hic_mols, &variant_contig_order);
    assess_breakpoints(&hic_links, &params, &variant_contig_order, &assembly);
}

fn good_assembly_loci(
    assembly: &Assembly,
    allele_fractions: &HashMap<i32, f32>,
    bad_alleles: &HashSet<i32>,
    
) -> ContigLoci {
    // returning a map from kmer id to contig id and position
    let mut variant_contig_order: HashMap<i32, ContigLocus> = HashMap::new();

    let mut contig_positions: HashMap<i32, Vec<(usize, i32, i32)>> = HashMap::new();
    for (kmer, (contig, num, _order, position)) in assembly.variants.iter() {
      
        if assembly.variants.contains_key(&Kmers::pair(*kmer)) {
            continue;
        } // we see both ref and alt in assembly, skip
        
        if let Some(fraction) = allele_fractions.get(&kmer.abs()) {
            if *fraction < MIN_ALLELE_FRACTION_HIC {
                continue;
            }
            
        } 
        if bad_alleles.contains(kmer) {
            continue;
        }

        if *num > 1 {
            continue;
        } // we see this kmer multiple times in the assembly, skip
        let positions = contig_positions.entry(*contig).or_insert(Vec::new());
        positions.push((*position, *kmer, *contig));
    }

    let mut loci: HashMap<i32, Vec<ContigLocus>> = HashMap::new();
    for (contig, positions) in contig_positions.iter() {
        let mut poses: Vec<(usize, i32, i32)> = Vec::new();
        for (pos, kmer, contig) in positions {
            poses.push((*pos, *kmer, *contig));
        }
        poses.sort();
        let contig_loci = loci.entry(*contig).or_insert(Vec::new());
        for (index, (position, kmer, contig)) in poses.iter().enumerate() {
            let reference;
            let alternate;
            match allele(*kmer) {
                Allele::Ref => {
                    reference = kmer.abs();
                    alternate = reference - 1;
                }
                Allele::Alt => {
                    alternate = kmer.abs();
                    reference = alternate + 1;
                }
            }
            let locus = ContigLocus {
                contig_id: *contig,
                index: index,
                position: *position,
                allele: allele(*kmer),
                reference: reference,
                alternate: alternate,
            };
            variant_contig_order.insert(*kmer, locus);
            variant_contig_order.insert(Kmers::pair(*kmer), locus);
            contig_loci.push(locus);
        }
    }

    ContigLoci {
        kmers: variant_contig_order,
        loci: loci,
    }
}

fn get_bad_alleles(hic_mols: &Mols) -> HashSet<i32> {
    let mut bad: HashSet<i32> = HashSet::new();
    for mol in hic_mols.get_molecules() {
        for i in 0..mol.len() {
            for j in (i + 1)..mol.len() {
                if Kmers::pair(mol[i].abs()) == mol[j].abs() {
                    bad.insert(mol[i].abs());
                    bad.insert(mol[j].abs());
                }
            }
        }
    }
    bad
}

fn get_allele_fractions(hic_mols: &Mols, kmers: &Kmers) -> HashMap<i32, f32> {
    let mut allele_fractions: HashMap<i32, f32> = HashMap::new();
    let mut allele_counts: HashMap<i32, [u32; 2]> = HashMap::new();
    for mol in hic_mols.get_molecules() {
        for var in mol {
            
            
            let canonical = var.abs().min(Kmers::pair(var.abs()));
            let count = allele_counts.entry(canonical).or_insert([0; 2]);
            if var.abs() % 2 == 0 {
                count[0] += 1;
            } else {
                count[1] += 1;
            }
        }
    }
    for (canonical, counts) in allele_counts {

        match kmers.kmer_type.get(&canonical).unwrap() {
            KmerType::PairedHet => (),
            KmerType::UnpairedHet => continue,
            KmerType::Homozygous => continue,
        }
        let fraction = (counts[0].min(counts[1]) as f32) / ((counts[0] + counts[1]) as f32);
        allele_fractions.insert(canonical, fraction);
        allele_fractions.insert(Kmers::pair(canonical), fraction);
    }
    allele_fractions
}

fn kmer_contig_position(kmer: i32, assembly: &Assembly, any: bool) -> Option<(i32, usize, usize)> {
    if let Some((contig_id, number_seen, order, position)) = assembly.variants.get(&kmer.abs()) {
        if any || *number_seen == 1 {
            return Some((*contig_id, *position, *order));
        }
    } else if let Some((contig_id, number_seen, order, position)) = assembly.variants.get(&Kmers::pair(kmer.abs())) {
        if any || *number_seen == 1 {
            return Some((*contig_id, *position, *order));
        }
    }
    None
}



fn gather_hic_links(
    hic_molecules: &Mols,
    variant_contig_order: &ContigLoci,
) -> HashMap<i32, Vec<Molecule>> {
    // returns map from contig id to list of HIC data structures
    let mut hic_mols: HashMap<i32, Vec<Molecule>> = HashMap::new();
    let mut long_hic_mols: HashMap<i32, Vec<Molecule>> = HashMap::new();
    let mut total = 0;
    for (contig, _) in variant_contig_order.loci.iter() {
        hic_mols.insert(*contig, Vec::new());
        long_hic_mols.insert(*contig, Vec::new());
    }

    let mut used_count = 0;
    let mut not_assembly = 0;
    let mut diff_contig = 0;
    for mol in hic_molecules.get_molecules() {
        let mut the_contig: Option<i32> = None;
        let mut loci: Vec<usize> = Vec::new();
        let mut alleles: Vec<Allele> = Vec::new();
        let mut used: HashSet<i32> = HashSet::new();
        let mut min: usize = std::usize::MAX;
        let mut max: usize = 0;

        for var in mol {
            if used.contains(&var.abs()) {
                used_count += 1;
                continue;
            }
            if let Some(locus) = variant_contig_order.kmers.get(&var.abs()) {
                if let Some(chrom) = the_contig {
                    min = min.min(locus.position);
                    max = max.max(locus.position);
                    if locus.contig_id == chrom {
                        loci.push(locus.index);
                        alleles.push(allele(*var));
                    } else {
                        diff_contig += 1;
                    }
                } else {
                    min = min.min(locus.position);
                    max = max.max(locus.position);
                    the_contig = Some(locus.contig_id);
                    loci.push(locus.index);
                    alleles.push(allele(*var));
                }
            } else {
                not_assembly += 1;
            }
            used.insert(var.abs());
        }
        if loci.len() > 1 {
            let contig_mols = hic_mols.entry(the_contig.unwrap()).or_insert(Vec::new());
           
            contig_mols.push(Molecule {
                loci: loci,
                alleles: alleles,
            });
            total += 1;
        }
    }
    eprintln!(
        "after culling we have {} hic molecules hitting >=2 distinct loci",
        total
    );
    eprintln!("why did we lose kmers? overlaps (same kmer twice) {}, no assembly locus {}, cross contig {}", used_count, not_assembly, diff_contig);

    hic_mols
}


fn allele(kmer: i32) -> Allele {
    match kmer.abs() % 2 == 0 {
        true => Allele::Alt,
        false => Allele::Ref,
    }
}

#[derive(Debug, Clone, Copy)]
enum Allele {
    Ref,
    Alt,
}

#[derive(Clone)]
struct Molecule {
    loci: Vec<usize>,
    alleles: Vec<Allele>,
}

fn check_remove(
    mol: &Molecule,
    mol_index: usize,
    current_mol_set: &mut HashSet<usize>,
    left: usize,
    right: usize,
) {
    let mut count = 0;
    for locus in mol.loci.iter() {
        if *locus >= left && *locus < right {
            count += 1;
        }
    }
    if count < 2 {
        current_mol_set.remove(&mol_index);
    }
}

fn check_add(
    mol: &Molecule,
    mol_index: usize,
    current_mol_set: &mut HashSet<usize>,
    left: usize,
    right: usize,
) {
    let mut count = 0;
    for locus in mol.loci.iter() {
        if *locus >= left && *locus < right {
            count += 1;
        }
    }
    if count < 2 {
        current_mol_set.insert(mol_index);
    }
}


/*
fn assess_breakpoints(
    hic_links: &HashMap<i32, Vec<Molecule>>,
    params: &Params,
    assembly: &Assembly,
) -> (
    HashMap<i32, Vec<(usize, usize)>>,
    HashMap<i32, Vec<(usize, usize)>>,
) {
    let mut chunks: HashMap<i32, Vec<(usize, usize)>> = HashMap::new(); // ranges for each contig
    let mut chunks_indices: HashMap<i32, Vec<(usize, usize)>> = HashMap::new(); // ranges for each contig



    //for (contig, hic) in hic_links.iter() {
    for contig in 1..assembly.contig_names.len() {
        let contig_name = &assembly.contig_names[contig];
        let contig = &(contig as i32);
        let empty: Vec<Molecule> = Vec::new();
        let hic = hic_links.get(contig).unwrap_or(&empty); //(&format!("contig {} {} has no hic links, contig size {}", contig, contig_name, assembly.contig_sizes.get(contig).unwrap()));
        let mut in_chunk = true;
        let contig_chunk = chunks.entry(*contig).or_insert(Vec::new());
        let contig_chunk_indices = chunks_indices.entry(*contig).or_insert(Vec::new());
        let mut current_chunk = (0, 0);
        let mut current_chunk_indices = (0, 0);
        let empty2: Vec<Option<bool>> = Vec::new();
        //let contig_phasing = phasing.phasing.get(contig).unwrap_or(&empty2); //expect("contig not in phasing?");

        let mut locus_hic_mols: HashMap<usize, Vec<usize>> = HashMap::new();
        for (read_index, hic_read) in hic.iter().enumerate() {
            for locus in hic_read.loci.iter() {
                let blah = locus_hic_mols.entry(*locus).or_insert(Vec::new());
                blah.push(read_index);
            }
        }
     
        let kmer_positions = assembly.contig_kmers.get(&contig).expect("contig not in assembly????");
        

        let mut current_hic_mol_set: HashSet<usize> = HashSet::new();
        for index in 0..params.break_window {
            let left = 0;
            let right = params.break_window;
            if let Some(hic_indexes) = locus_hic_mols.get(&index) {
                for hic_index in hic_indexes {
                    check_add(
                        &hic[*hic_index],
                        *hic_index,
                        &mut current_hic_mol_set,
                        left,
                        right,
                    );
                }
            }
        }
        for (middle_index, (position, kmer)) in kmer_positions.iter().enumerate() {
            let mut left = middle_index - params.break_window;
            let right = middle_index + params.break_window;
            if middle_index > params.break_window {
                if let Some(hic_indexes) = locus_hic_mols.get(&(left - 1)) {
                    for hic_index in hic_indexes {
                        check_remove(
                            &hic[*hic_index],
                            *hic_index,
                            &mut current_hic_mol_set,
                            left,
                            right,
                        );
                    }
                }
            } else {
                left = 0;
            }
            if let Some(hic_indexes) = locus_hic_mols.get(&(right - 1)) {
                for hic_index in hic_indexes {
                    check_add(
                        &hic[*hic_index],
                        *hic_index,
                        &mut current_hic_mol_set,
                        left,
                        right,
                    );
                }
            }

            let mut hic_links = 0;
            if middle_index > 250 && middle_index < kmer_positions.len() - 250 {
                for hic_moldex in current_hic_mol_set.iter() {
                    let hicmol = &hic[*hic_moldex];
                    let mut molcounts: [u16; 4] = [0; 4];
                    for index1 in 0..hicmol.loci.len() {
                        let locus1 = hicmol.loci[index1];
                        for index2 in (index1 + 1)..hicmol.loci.len() {
                            let locus2 = hicmol.loci[index2];
                            if locus1 < middle_index && locus1 >= left && locus2 >= middle_index && locus2 < right {
                                hic_links += 1;
                                
                            }
                        }
                    }
                }

            
                if hic_links > 10 && in_chunk {
                    current_chunk = (current_chunk.0, *position);
                    current_chunk_indices = (current_chunk_indices.0, middle_index);
                } else if hic_links <= 10 && in_chunk {
                    in_chunk = false;
                    if current_chunk.1 > current_chunk.0 {
                        contig_chunk.push(current_chunk);
                        contig_chunk_indices.push(current_chunk_indices);
                        eprintln!(
                            "adding chunk for contig, total hic links {}, contig {}, chunk {:?}",
                            hic_links, contig_name, current_chunk
                        );
                    }
                    current_chunk = (position + 1, position + 1);
                    current_chunk_indices = (middle_index + 1, middle_index + 1);
                } else if hic_links > 10 && !in_chunk {
                    in_chunk = true;
                    //current_chunk = (locus.position, locus.position);
                }
            } else if in_chunk {
                current_chunk = (current_chunk.0, *position);
                current_chunk_indices = (current_chunk_indices.0, middle_index);
            }
        }
        if in_chunk || contig_chunk.len() == 0 {
            current_chunk = (current_chunk.0, *assembly.contig_sizes.get(contig).unwrap());
            contig_chunk.push(current_chunk);
            
            current_chunk_indices = (
                current_chunk_indices.0,
                kmer_positions.len(),
            );
            contig_chunk_indices.push(current_chunk_indices);
            eprintln!(
                "adding chunk at finish for contig {}, chunk {:?}, and indices {:?}",
                contig_name, current_chunk, current_chunk_indices
            );
        }
        if contig_chunk.len() > 1 {
            eprintln!(
                "contig {} with size {} is split into {} chunks",
                contig,
                *assembly.contig_sizes.get(contig).unwrap(),
                contig_chunk.len()
            );
            for (start, end) in contig_chunk.iter() {
                eprintln!("\t{} - {}", start, end);
            }
            for (startdex, enddex) in contig_chunk_indices.iter() {
                eprintln!("\t{} - {}", startdex, enddex);
            }
        } else {
            eprintln!("contig {} has no breaks", contig);
            for (start, end) in contig_chunk.iter() {
                eprintln!("\t{} - {}", start, end);
            }
            for (startdex, enddex) in contig_chunk_indices.iter() {
                eprintln!("\t{} - {}", startdex, enddex);
            }
        }
    }
    let reader = fasta::Reader::from_file(Path::new(&params.assembly_fasta.to_string()))
        .expect("fasta not found");
    let mut writer = fasta::Writer::to_file(Path::new(&format!("{}/breaks.fa", params.output)))
        .expect("cannot open fasta writer");
    for record in reader.records() {
        let record = record.unwrap();
        let contig_name = record.id().to_string();
        let contig_id = assembly.contig_ids.get(&contig_name).unwrap();

        if !chunks.contains_key(contig_id) {
            eprintln!("contig has no chunks??? {}", contig_id);
            let size = assembly.contig_sizes.get(contig_id).unwrap();
            let range = chunks.entry(*contig_id).or_insert(Vec::new());
            range.push((0, *size));
        }
        let ranges = chunks.get(contig_id).unwrap();

        for (index, (start, stop)) in ranges.iter().enumerate() {
            let mut new_contig_name = contig_name.to_string();
            if ranges.len() > 0 {
                let list = vec![
                    new_contig_name,
                    (index + 1).to_string(),
                    start.to_string(),
                    stop.to_string(),
                ];
                new_contig_name = list.join("_");
            }
            let seq: TextSlice = &record.seq()[*start..*stop];
            let record = Record::with_attrs(&new_contig_name, None, &seq);
            writer
                .write_record(&record)
                .expect("could not write record");
        }
    }
    (chunks, chunks_indices)
}
*/


fn assess_breakpoints(
    hic_links: &HashMap<i32, Vec<Molecule>>,
    params: &Params,
    contig_loci: &ContigLoci,
    assembly: &Assembly,
) -> (
    HashMap<i32, Vec<(usize, usize)>>,
    HashMap<i32, Vec<(usize, usize)>>,
) {
    let mut chunks: HashMap<i32, Vec<(usize, usize)>> = HashMap::new(); // ranges for each contig
    let mut chunks_indices: HashMap<i32, Vec<(usize, usize)>> = HashMap::new(); // ranges for each contig

    //for (contig, hic) in hic_links.iter() {
    for contig in 1..assembly.contig_names.len() {
        let contig_name = &assembly.contig_names[contig];
        let contig = &(contig as i32);
        let empty: Vec<Molecule> = Vec::new();
        let hic = hic_links.get(contig).unwrap_or(&empty); //(&format!("contig {} {} has no hic links, contig size {}", contig, contig_name, assembly.contig_sizes.get(contig).unwrap()));
        let mut in_chunk = true;
        let contig_chunk = chunks.entry(*contig).or_insert(Vec::new());
        let contig_chunk_indices = chunks_indices.entry(*contig).or_insert(Vec::new());
        let mut current_chunk = (0, 0);
        let mut current_chunk_indices = (0, 0);
        let empty2: Vec<Option<bool>> = Vec::new();

        let mut locus_hic_mols: HashMap<usize, Vec<usize>> = HashMap::new();
        for (read_index, hic_read) in hic.iter().enumerate() {
            for locus in hic_read.loci.iter() {
                let blah = locus_hic_mols.entry(*locus).or_insert(Vec::new());
                blah.push(read_index);
            }
        }
        let empty3: Vec<ContigLocus> = Vec::new();
        let loci = contig_loci.loci.get(contig).unwrap_or(&empty3);

        let mut current_hic_mol_set: HashSet<usize> = HashSet::new();
        for index in 0..params.break_window {
            let left = 0;
            let right = params.break_window;
            if let Some(hic_indexes) = locus_hic_mols.get(&index) {
                for hic_index in hic_indexes {
                    check_add(
                        &hic[*hic_index],
                        *hic_index,
                        &mut current_hic_mol_set,
                        left,
                        right,
                    );
                }
            }
        }
        for (mid, locus) in loci.iter().enumerate() {
            let mut left = mid - params.break_window;
            let right = mid + params.break_window;
            if mid > params.break_window {
                if let Some(hic_indexes) = locus_hic_mols.get(&(left - 1)) {
                    for hic_index in hic_indexes {
                        check_remove(
                            &hic[*hic_index],
                            *hic_index,
                            &mut current_hic_mol_set,
                            left,
                            right,
                        );
                    }
                }
            } else {
                left = 0;
            }
            if let Some(hic_indexes) = locus_hic_mols.get(&(right - 1)) {
                for hic_index in hic_indexes {
                    check_add(
                        &hic[*hic_index],
                        *hic_index,
                        &mut current_hic_mol_set,
                        left,
                        right,
                    );
                }
            }
            let mut hic_links = 0;

            for hic_moldex in current_hic_mol_set.iter() {
                let hicmol = &hic[*hic_moldex];
                let mut molcounts: [u16; 4] = [0; 4];
                for index1 in 0..hicmol.loci.len() {
                    for index2 in (index1 + 1)..hicmol.loci.len() {
                        let locus1 = hicmol.loci[index1];
                        let locus2 = hicmol.loci[index2];
                        if locus1 < mid && locus1 >= left && locus2 >= mid && locus2 < right {
                            hic_links += 1;
                            
                        }
                    }
                }
            }
            //eprintln!(
            //    "contig {}, mid {} position {}, break_counts {:?}",
            //   contig, locus.position, mid, counts
            //);
            if mid > 250 && mid < loci.len() - 250 {
                if hic_links > 10 && in_chunk {
                    current_chunk = (current_chunk.0, locus.position);
                    current_chunk_indices = (current_chunk_indices.0, mid);
                } else if hic_links <= 10 && in_chunk {
                    in_chunk = false;
                    if current_chunk.1 > current_chunk.0 {
                        contig_chunk.push(current_chunk);
                        contig_chunk_indices.push(current_chunk_indices);
                        eprintln!(
                            "adding chunk for contig, hic_links {} contig {}, chunk {:?}",
                            hic_links,  contig_name, current_chunk
                        );
                    }
                    current_chunk = (locus.position + 1, locus.position + 1);
                    current_chunk_indices = (mid + 1, mid + 1);
                } else if hic_links > 10 && !in_chunk {
                    in_chunk = true;
                    //current_chunk = (locus.position, locus.position);
                }
            } else if in_chunk {
                current_chunk = (current_chunk.0, locus.position);
                current_chunk_indices = (current_chunk_indices.0, mid);
            }
        }
        if in_chunk || contig_chunk.len() == 0 {
            current_chunk = (current_chunk.0, *assembly.contig_sizes.get(contig).unwrap());
            contig_chunk.push(current_chunk);
            let empty: Vec<ContigLocus> = Vec::new();
            current_chunk_indices = (
                current_chunk_indices.0,
                contig_loci.loci.get(contig).unwrap_or(&empty).len(),
            );
            contig_chunk_indices.push(current_chunk_indices);
            eprintln!(
                "adding chunk at finish for contig {}, chunk {:?}, and indices {:?}",
                contig_name, current_chunk, current_chunk_indices
            );
        }
        if contig_chunk.len() > 1 {
            eprintln!(
                "contig {} with size {} is split into {} chunks",
                contig,
                *assembly.contig_sizes.get(contig).unwrap(),
                contig_chunk.len()
            );
            for (start, end) in contig_chunk.iter() {
                eprintln!("\t{} - {}", start, end);
            }
            for (startdex, enddex) in contig_chunk_indices.iter() {
                eprintln!("\t{} - {}", startdex, enddex);
            }
        } else {
            eprintln!("contig {} has no breaks", contig);
            for (start, end) in contig_chunk.iter() {
                eprintln!("\t{} - {}", start, end);
            }
            for (startdex, enddex) in contig_chunk_indices.iter() {
                eprintln!("\t{} - {}", startdex, enddex);
            }
        }
    }
    let reader = fasta::Reader::from_file(Path::new(&params.assembly_fasta.to_string()))
        .expect("fasta not found");
    let mut writer = fasta::Writer::to_file(Path::new(&format!("{}/breaks.fa", params.output)))
        .expect("cannot open fasta writer");
    for record in reader.records() {
        let record = record.unwrap();
        let contig_name = record.id().to_string();
        let contig_id = assembly.contig_ids.get(&contig_name).unwrap();

        if !chunks.contains_key(contig_id) {
            eprintln!("contig has no chunks??? {}", contig_id);
            let size = assembly.contig_sizes.get(contig_id).unwrap();
            let range = chunks.entry(*contig_id).or_insert(Vec::new());
            range.push((0, *size));
        }
        let ranges = chunks.get(contig_id).unwrap();

        for (index, (start, stop)) in ranges.iter().enumerate() {
            let mut new_contig_name = contig_name.to_string();
            if ranges.len() > 0 {
                let list = vec![
                    new_contig_name,
                    (index + 1).to_string(),
                    start.to_string(),
                    stop.to_string(),
                ];
                new_contig_name = list.join("_");
            }
            let seq: TextSlice = &record.seq()[*start..*stop];
            let record = Record::with_attrs(&new_contig_name, None, &seq);
            writer
                .write_record(&record)
                .expect("could not write record");
        }
    }
    (chunks, chunks_indices)
}



#[derive(Clone)]
struct Params {
    het_kmers: String,
    hic_mols: Vec<String>,
    min_contig_length: usize,
    output: String,
    assembly_kmers: String,
    assembly_fasta: String,
    threads: usize,
    break_window: usize,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();

    let het_kmers = params.value_of("het_kmers").unwrap();
    let output = params.value_of("output").unwrap();

    let mut hic_mols: Vec<String> = Vec::new();
    match params.value_of("hic_mols") {
        Some(hic_fofn) => {
            let f = File::open(hic_fofn).expect("Unable to open hic fofn");
            let f = BufReader::new(f);

            for line in f.lines() {
                let line = line.expect("Unable to read txg fofn line");
                hic_mols.push(line.to_string());
            }
        },
        None => (),
    }

    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<usize>().unwrap();

    let assembly_kmers = params.value_of("assembly_kmers").unwrap();
    let assembly_fasta = params.value_of("assembly_fasta").unwrap();

    let break_window = params.value_of("break_window").unwrap_or("500");
    let break_window = break_window.to_string().parse::<usize>().unwrap();
    eprintln!("break window {}", break_window);

    let min_contig_length = params.value_of("min_contig_length").unwrap_or("100000");
    let min_contig_length = min_contig_length.to_string().parse::<usize>().unwrap();

    Params {
        het_kmers: het_kmers.to_string(),
        output: output.to_string(),
        hic_mols: hic_mols,
        assembly_kmers: assembly_kmers.to_string(),
        assembly_fasta: assembly_fasta.to_string(),
        threads: threads,
        break_window: break_window,
        min_contig_length: min_contig_length,
    }
}
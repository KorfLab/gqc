package main

import (
	"bufio"
	"compress/gzip"
	"encoding/gob"
	// "flag"
	"fmt"
	"io"
	"math"
	"net/http"
	"net/url"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/PuerkitoBio/goquery"
)

// -----------------------------------------------------------------------------
// UTILITIES & CONSTANTS
// -----------------------------------------------------------------------------

func getReader(filename string) (io.ReadCloser, *bufio.Scanner) {
	file, err := os.Open(filename)
	if err != nil {
		panic(fmt.Sprintf("Cannot open %s: %v", filename, err))
	}
	if strings.HasSuffix(filename, ".gz") {
		gz, err := gzip.NewReader(file)
		if err != nil {
			panic(err)
		}
		// Wrap in a custom ReadCloser to close both gz and file
		rc := struct {
			io.Reader
			io.Closer
		}{gz, file}
		return rc, bufio.NewScanner(gz)
	}
	return file, bufio.NewScanner(file)
}

func prob2score(p float64) float64 {
	if p == 0.0 {
		return -100.0
	}
	return math.Log2(p / 0.25)
}

func anti(seq string) string {
	b := []byte(seq)
	n := len(b)
	res := make([]byte, n)
	for i := 0; i < n; i++ {
		c := b[n-1-i]
		switch c {
		case 'A', 'a': res[i] = 'T'
		case 'C', 'c': res[i] = 'G'
		case 'G', 'g': res[i] = 'C'
		case 'T', 't': res[i] = 'A'
		case 'R', 'r': res[i] = 'Y'
		case 'Y', 'y': res[i] = 'R'
		case 'K', 'k': res[i] = 'M'
		case 'M', 'm': res[i] = 'K'
		case 'W', 'w': res[i] = 'W'
		case 'S', 's': res[i] = 'S'
		case 'B', 'b': res[i] = 'V'
		case 'V', 'v': res[i] = 'B'
		case 'D', 'd': res[i] = 'H'
		case 'H', 'h': res[i] = 'D'
		default: res[i] = c
		}
	}
	return string(res)
}

// -----------------------------------------------------------------------------
// DATA STRUCTURES & CACHING
// -----------------------------------------------------------------------------

type Feature struct {
	Seqid  string
	Typ    string
	Beg    int
	End    int
	Strand byte
	Fid    string
	Pid    string
}

type Transcript struct {
	TxFeat   Feature
	Exons    []Feature
	Introns  []Feature
	Cdss     []Feature
	Utr5s    []Feature
	Utr3s    []Feature
	IsCoding bool
}

func (t *Transcript) Finalize() {
	t.IsCoding = len(t.Cdss) > 0
}

type GeneGroup struct {
	Gene Feature
	Txs  []Transcript
}

type GenomeData struct {
	Sequences map[string]string
	Genes     []GeneGroup
}

func (gd *GenomeData) LoadFasta(path string, v bool) {
	if v {
		fmt.Printf("  [Verbose] Parsing FASTA: %s\n", path)
	}
	rc, scanner := getReader(path)
	defer rc.Close()

	if gd.Sequences == nil {
		gd.Sequences = make(map[string]string)
	}

	var name string
	var seqBuilder strings.Builder

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if strings.HasPrefix(line, ">") {
			if name != "" {
				gd.Sequences[name] = seqBuilder.String()
				seqBuilder.Reset()
			}
			parts := strings.Fields(line[1:])
			name = strings.TrimPrefix(parts[0], "chr")
		} else {
			seqBuilder.WriteString(line)
		}
	}
	if name != "" {
		gd.Sequences[name] = seqBuilder.String()
	}
}

func (gd *GenomeData) LoadGff3(path string, v bool) {
	if v {
		fmt.Printf("  [Verbose] Parsing GFF3: %s\n", path)
	}
	rc, scanner := getReader(path)
	defer rc.Close()

	featuresByPid := make(map[string][]Feature)
	var genesRaw []Feature

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") || len(line) == 0 {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) != 9 {
			continue
		}

		seqid := strings.TrimPrefix(fields[0], "chr")
		typ := fields[2]
		beg, _ := strconv.Atoi(fields[3])
		end, _ := strconv.Atoi(fields[4])
		strand := fields[6][0]

		var fid string
		pids := []string{""}

		attrs := strings.Split(strings.TrimRight(fields[8], ";"), ";")
		for _, attr := range attrs {
			parts := strings.SplitN(attr, "=", 2)
			if len(parts) == 2 {
				if parts[0] == "ID" {
					fid = parts[1]
				} else if parts[0] == "Parent" {
					pids = strings.Split(parts[1], ",")
				}
			}
		}

		for _, pid := range pids {
			feat := Feature{Seqid: seqid, Typ: typ, Beg: beg, End: end, Strand: strand, Fid: fid, Pid: pid}
			if typ == "gene" {
				genesRaw = append(genesRaw, feat)
			} else {
				featuresByPid[pid] = append(featuresByPid[pid], feat)
			}
		}
	}

	for _, geneFeat := range genesRaw {
		var txs []Transcript
		for _, txFeat := range featuresByPid[geneFeat.Fid] {
			tx := Transcript{TxFeat: txFeat}
			for _, child := range featuresByPid[txFeat.Fid] {
				switch child.Typ {
				case "exon":
					tx.Exons = append(tx.Exons, child)
				case "intron":
					tx.Introns = append(tx.Introns, child)
				case "CDS":
					tx.Cdss = append(tx.Cdss, child)
				case "five_prime_UTR":
					tx.Utr5s = append(tx.Utr5s, child)
				case "three_prime_UTR":
					tx.Utr3s = append(tx.Utr3s, child)
				}
			}
			tx.Finalize()
			txs = append(txs, tx)
		}
		gd.Genes = append(gd.Genes, GeneGroup{Gene: geneFeat, Txs: txs})
	}
}

func buildOrLoadCache(dbPath string, fastaPath string, gff3Path string, force bool, v bool) *GenomeData {
	cacheFile := fmt.Sprintf("%s.bin", dbPath)
	if !force {
		if _, err := os.Stat(cacheFile); err == nil {
			if v {
				fmt.Printf("  [Verbose] Loading cached database from %s\n", cacheFile)
			}
			file, _ := os.Open(cacheFile)
			defer file.Close()
			reader := bufio.NewReader(file)
			decoder := gob.NewDecoder(reader)
			var gd GenomeData
			decoder.Decode(&gd)
			return &gd
		}
	}

	gd := &GenomeData{}
	gd.LoadFasta(fastaPath, v)
	gd.LoadGff3(gff3Path, v)

	if v {
		fmt.Printf("  [Verbose] Serializing database to %s\n", cacheFile)
	}
	file, _ := os.Create(cacheFile)
	defer file.Close()
	writer := bufio.NewWriter(file)
	encoder := gob.NewEncoder(writer)
	encoder.Encode(gd)
	writer.Flush()

	return gd
}

// -----------------------------------------------------------------------------
// MARKOV MODEL
// -----------------------------------------------------------------------------

type Model struct {
	K            int
	Probs        map[string]map[byte]float64
	LogOddsCache map[string]float64
}

func NewModel() *Model {
	return &Model{
		Probs:        make(map[string]map[byte]float64),
		LogOddsCache: make(map[string]float64),
	}
}

func (m *Model) Create(seqs []string, k int) {
	m.K = k
	counts := make(map[string][]uint32)

	ntToIdx := func(b byte) int {
		switch b {
		case 'A', 'a': return 0
		case 'C', 'c': return 1
		case 'G', 'g': return 2
		case 'T', 't': return 3
		}
		return -1
	}

	for _, seq := range seqs {
		if len(seq) <= k {
			continue
		}
		for i := k; i < len(seq); i++ {
			idx := ntToIdx(seq[i])
			if idx != -1 {
				kmer := seq[i-k : i]
				if counts[kmer] == nil {
					counts[kmer] = []uint32{1, 1, 1, 1} // Pseudocounts
				}
				counts[kmer][idx]++
			}
		}
	}

	idxToNt := []byte{'A', 'C', 'G', 'T'}
	for kmer, cnts := range counts {
		var total uint32
		for _, v := range cnts {
			total += v
		}
		m.Probs[kmer] = make(map[byte]float64)
		for i := 0; i < 4; i++ {
			p := float64(cnts[i]) / float64(total)
			m.Probs[kmer][idxToNt[i]] = p
			cacheKey := kmer + string(idxToNt[i])
			m.LogOddsCache[cacheKey] = prob2score(p)
		}
	}
}

func (m *Model) Score(seq string) float64 {
	if len(seq) <= m.K {
		return 0.0
	}
	var total float64
	for i := 0; i <= len(seq)-m.K-1; i++ {
		window := seq[i : i+m.K+1]
		if val, exists := m.LogOddsCache[window]; exists {
			total += val
		} else {
			total += -2.0
		}
	}
	return total
}

func (m *Model) ExportFile(path string) {
	file, _ := os.Create(path)
	defer file.Close()
	writer := bufio.NewWriter(file)

	fmt.Fprintf(writer, "%% MM %s %d\n", path, len(m.Probs)*4)

	var keys []string
	for k := range m.Probs {
		keys = append(keys, k)
	}
	sort.Strings(keys)

	for _, kmer := range keys {
		nts := m.Probs[kmer]
		for _, nt := range []byte{'A', 'C', 'G', 'T'} {
			if p, exists := nts[nt]; exists {
				fmt.Fprintf(writer, "%s%c %.6f\n", kmer, nt, p)
			}
		}
		fmt.Fprintln(writer)
	}
	writer.Flush()
}

func (m *Model) ImportFile(path string) {
	rc, scanner := getReader(path)
	defer rc.Close()

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "%") || strings.TrimSpace(line) == "" {
			continue
		}
		parts := strings.Fields(line)
		if len(parts) != 2 || len(parts[0]) < 2 {
			continue
		}

		key := parts[0]
		kmer := key[:len(key)-1]
		nt := key[len(key)-1]
		val, _ := strconv.ParseFloat(parts[1], 64)

		if m.K == 0 {
			m.K = len(kmer)
		}

		if m.Probs[kmer] == nil {
			m.Probs[kmer] = make(map[byte]float64)
		}
		m.Probs[kmer][nt] = val
		m.LogOddsCache[key] = prob2score(val)
	}
}

// -----------------------------------------------------------------------------
// ASYNC BULK PIPELINE (The Auto Crawler)
// -----------------------------------------------------------------------------

type GenomeJob struct {
	Prefix   string
	FastaURL string
	GffURL   string
}

func crawlDatabase(baseURLStr string) ([]GenomeJob, error) {
	fmt.Printf("Crawling %s ...\n", baseURLStr)
	baseURL, err := url.Parse(baseURLStr)
	if err != nil {
		return nil, err
	}

	// 1. Create a custom request to spoof the User-Agent
	req, err := http.NewRequest("GET", baseURLStr, nil)
	if err != nil {
		return nil, err
	}
	// NCBI requires custom user agents to avoid blocking generic bots
	req.Header.Set("User-Agent", "GeneQC-Pipeline/1.0 (Research Tool)")

	client := &http.Client{Timeout: 30 * time.Second}
	res, err := client.Do(req)
	if err != nil {
		return nil, err
	}
	defer res.Body.Close()

	// 2. Catch 403 Forbidden or 404 Not Found errors loudly!
	if res.StatusCode != 200 {
		return nil, fmt.Errorf("server rejected the request (HTTP %d %s)", res.StatusCode, res.Status)
	}

	doc, err := goquery.NewDocumentFromReader(res.Body)
	if err != nil {
		return nil, err
	}

	var fastas []string
	var gffs []string

	doc.Find("a").Each(func(i int, s *goquery.Selection) {
		href, exists := s.Attr("href")
		if exists {
			// Safely parse the href rather than assuming it's a raw path
			parsedHref, parseErr := url.Parse(href)
			if parseErr == nil {
				fullURL := baseURL.ResolveReference(parsedHref).String()
				if strings.HasSuffix(href, ".fna.gz") || strings.HasSuffix(href, ".fa.gz") {
					fastas = append(fastas, fullURL)
				} else if strings.HasSuffix(href, ".gff.gz") || strings.HasSuffix(href, ".gff3.gz") {
					gffs = append(gffs, fullURL)
				}
			}
		}
	})

	var jobs []GenomeJob
	for _, fastaURL := range fastas {
		parts := strings.Split(fastaURL, "/")
		fName := parts[len(parts)-1]
		prefix := strings.Replace(fName, ".fna.gz", "", 1)
		prefix = strings.Replace(prefix, ".fa.gz", "", 1)

		for _, gffURL := range gffs {
			if strings.Contains(gffURL, prefix) {
				jobs = append(jobs, GenomeJob{Prefix: prefix, FastaURL: fastaURL, GffURL: gffURL})
				break
			}
		}
	}
	return jobs, nil
}

func downloadFile(fileURL string, dest string) error {
	out, err := os.Create(dest)
	if err != nil {
		return err
	}
	defer out.Close()

	// Use custom User-Agent for downloads as well
	req, err := http.NewRequest("GET", fileURL, nil)
	if err != nil {
		return err
	}
	req.Header.Set("User-Agent", "GeneQC-Pipeline/1.0 (Research Tool)")

	client := &http.Client{Timeout: 30 * time.Minute}
	resp, err := client.Do(req)
	if err != nil {
		return err
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return fmt.Errorf("bad status during download: %s", resp.Status)
	}

	_, err = io.Copy(out, resp.Body)
	return err
}

type BadGeneRes struct {
	Gid  string
	Tid  string
	Diff float64
}

func processGenomeWorker(job GenomeJob, k int, v bool, wg *sync.WaitGroup, sem chan struct{}) {
	defer wg.Done()
	defer func() { <-sem }() // Release concurrency token when done

	start := time.Now()
	fmt.Printf(">> Starting Genome: %s\n", job.Prefix)
	os.MkdirAll("bulk_out", os.ModePerm)

	fastaPath := fmt.Sprintf("%s_temp.fa.gz", job.Prefix)
	gffPath := fmt.Sprintf("%s_temp.gff.gz", job.Prefix)
	dbPath := fmt.Sprintf("%s_db", job.Prefix)

	// Clean up temp files automatically at the end of the function
	defer os.Remove(fastaPath)
	defer os.Remove(gffPath)
	defer os.Remove(fmt.Sprintf("%s.bin", dbPath))

	downloadFile(job.FastaURL, fastaPath)
	downloadFile(job.GffURL, gffPath)

	gd := buildOrLoadCache(dbPath, fastaPath, gffPath, true, v)

	mmbase := NewModel()
	var allSeqs []string
	for _, seq := range gd.Sequences {
		allSeqs = append(allSeqs, seq)
	}
	mmbase.Create(allSeqs, k)

	models := make(map[string]*Model)
	configs := []struct {
		name     string
		accessor func(*Transcript) []Feature
	}{
		{"exon", func(tx *Transcript) []Feature { return tx.Exons }},
		{"intron", func(tx *Transcript) []Feature { return tx.Introns }},
		{"utr5", func(tx *Transcript) []Feature { return tx.Utr5s }},
		{"utr3", func(tx *Transcript) []Feature { return tx.Utr3s }},
	}

	for _, cfg := range configs {
		var featSeqs []string
		for _, gGroup := range gd.Genes {
			for _, tx := range gGroup.Txs {
				for _, feat := range cfg.accessor(&tx) {
					seq, exists := gd.Sequences[feat.Seqid]
					if exists && feat.Beg > 0 && feat.End <= len(seq) {
						fSeq := seq[feat.Beg-1 : feat.End]
						if feat.Strand == '-' {
							fSeq = anti(fSeq)
						}
						featSeqs = append(featSeqs, fSeq)
					}
				}
			}
		}
		if len(featSeqs) > 0 {
			mm := NewModel()
			mm.Create(featSeqs, k)
			models[cfg.name] = mm
		}
	}

	var badgenes []BadGeneRes
	for _, gGroup := range gd.Genes {
		for _, tx := range gGroup.Txs {
			if !tx.IsCoding {
				continue
			}

			var tModelSc, tBaseSc float64
			var tModelPos, tBasePos int

			for _, cfg := range configs {
				if model, exists := models[cfg.name]; exists {
					for _, feat := range cfg.accessor(&tx) {
						seq, ok := gd.Sequences[feat.Seqid]
						if ok && feat.Beg > 0 && feat.End <= len(seq) {
							fSeq := seq[feat.Beg-1 : feat.End]
							if feat.Strand == '-' {
								fSeq = anti(fSeq)
							}

							tModelSc += model.Score(fSeq)
							pos1 := len(fSeq) - model.K
							if pos1 > 0 {
								tModelPos += pos1
							}

							tBaseSc += mmbase.Score(fSeq)
							pos2 := len(fSeq) - mmbase.K
							if pos2 > 0 {
								tBasePos += pos2
							}
						}
					}
				}
			}

			if tModelPos > 0 && tBasePos > 0 {
				diff := (tBaseSc / float64(tBasePos)) - (tModelSc / float64(tModelPos))
				if diff > 0.0 {
					badgenes = append(badgenes, BadGeneRes{Gid: gGroup.Gene.Fid, Tid: tx.TxFeat.Fid, Diff: diff})
				}
			}
		}
	}

	sort.Slice(badgenes, func(i, j int) bool { return badgenes[i].Diff > badgenes[j].Diff })

	outFile := fmt.Sprintf("bulk_out/%s_badgenes.txt", job.Prefix)
	f, _ := os.Create(outFile)
	defer f.Close()
	writer := bufio.NewWriter(f)
	fmt.Fprintln(writer, "GeneID\tTranscriptID\tScoreDiff")
	for _, bg := range badgenes {
		fmt.Fprintf(writer, "%s\t%s\t%.4f\n", bg.Gid, bg.Tid, bg.Diff)
	}
	writer.Flush()

	fmt.Printf("<< Analysis complete. Found %d anomalous transcripts. Results written to %s\n", len(badgenes), outFile)
	fmt.Printf("<< Finished Genome: %s in %v\n", job.Prefix, time.Since(start))
}

// -----------------------------------------------------------------------------
// CLI ENTRY POINT
// -----------------------------------------------------------------------------

func main() {
	// Parse global flags manually to allow git-style subcommands
	verbose := false
	args := os.Args[1:]

	if len(args) > 0 && (args[0] == "-v" || args[0] == "--verbose") {
		verbose = true
		args = args[1:]
	}

	if len(args) == 0 {
		fmt.Println("Usage: gene_qc [-v] <command> [args...]")
		fmt.Println("Commands: createdb, train, findbgs, auto")
		os.Exit(1)
	}

	command := args[0]
	cmdArgs := args[1:]

	switch command {
	case "createdb":
		if len(cmdArgs) < 3 {
			fmt.Println("Usage: gene_qc createdb <db_name> <fasta> <gff3>")
			return
		}
		fmt.Println("Building database...")
		start := time.Now()
		buildOrLoadCache(cmdArgs[0], cmdArgs[1], cmdArgs[2], true, verbose)
		fmt.Printf("Database saved to %s.bin in %v\n", cmdArgs[0], time.Since(start))

	case "train":
		if len(cmdArgs) < 1 {
			fmt.Println("Usage: gene_qc train <db_name> [-k <kmer_size>]")
			return
		}
		dbName := cmdArgs[0]
		k := 3
		if len(cmdArgs) == 3 && cmdArgs[1] == "-k" {
			k, _ = strconv.Atoi(cmdArgs[2])
		}

		fmt.Println("Loading database...")
		gd := buildOrLoadCache(dbName, "", "", false, verbose)
		os.MkdirAll("models", os.ModePerm)

		fmt.Println("Training baseline model...")
		mmbase := NewModel()
		var allSeqs []string
		for _, seq := range gd.Sequences {
			allSeqs = append(allSeqs, seq)
		}
		mmbase.Create(allSeqs, k)
		mmbase.ExportFile("models/baseline.mm")

		configs := []struct {
			name     string
			file     string
			accessor func(*Transcript) []Feature
		}{
			{"exon", "exon.mm", func(tx *Transcript) []Feature { return tx.Exons }},
			{"intron", "intron.mm", func(tx *Transcript) []Feature { return tx.Introns }},
			{"5' UTR", "5utr.mm", func(tx *Transcript) []Feature { return tx.Utr5s }},
			{"3' UTR", "3utr.mm", func(tx *Transcript) []Feature { return tx.Utr3s }},
		}

		for _, cfg := range configs {
			fmt.Printf("Training %s model...\n", cfg.name)
			var featSeqs []string
			for _, gGroup := range gd.Genes {
				for _, tx := range gGroup.Txs {
					for _, feat := range cfg.accessor(&tx) {
						seq, ok := gd.Sequences[feat.Seqid]
						if ok && feat.Beg > 0 && feat.End <= len(seq) {
							fSeq := seq[feat.Beg-1 : feat.End]
							if feat.Strand == '-' {
								fSeq = anti(fSeq)
							}
							featSeqs = append(featSeqs, fSeq)
						}
					}
				}
			}
			if len(featSeqs) > 0 {
				mm := NewModel()
				mm.Create(featSeqs, k)
				mm.ExportFile(filepath.Join("models", cfg.file))
			}
		}
		fmt.Println("Training complete! Models saved in 'models/' directory.")

	case "findbgs":
		if len(cmdArgs) < 1 {
			fmt.Println("Usage: gene_qc findbgs <db_name>")
			return
		}
		dbName := cmdArgs[0]
		fmt.Println("Loading database...")
		gd := buildOrLoadCache(dbName, "", "", false, verbose)

		fmt.Println("Loading models...")
		models := make(map[string]*Model)
		modelFiles := map[string]string{
			"base":   "baseline.mm",
			"exon":   "exon.mm",
			"intron": "intron.mm",
			"utr5":   "5utr.mm",
			"utr3":   "3utr.mm",
		}

		for key, file := range modelFiles {
			path := filepath.Join("models", file)
			if _, err := os.Stat(path); err == nil {
				m := NewModel()
				m.ImportFile(path)
				models[key] = m
			}
		}

		baseModel := models["base"]
		var badgenes []BadGeneRes

		fmt.Println("Scoring genes...")
		configs := []struct {
			name     string
			accessor func(*Transcript) []Feature
		}{
			{"exon", func(tx *Transcript) []Feature { return tx.Exons }},
			{"intron", func(tx *Transcript) []Feature { return tx.Introns }},
			{"utr5", func(tx *Transcript) []Feature { return tx.Utr5s }},
			{"utr3", func(tx *Transcript) []Feature { return tx.Utr3s }},
		}

		for _, gGroup := range gd.Genes {
			for _, tx := range gGroup.Txs {
				if !tx.IsCoding {
					continue
				}
				var tModelSc, tBaseSc float64
				var tModelPos, tBasePos int

				for _, cfg := range configs {
					if model, exists := models[cfg.name]; exists {
						for _, feat := range cfg.accessor(&tx) {
							seq, ok := gd.Sequences[feat.Seqid]
							if ok && feat.Beg > 0 && feat.End <= len(seq) {
								fSeq := seq[feat.Beg-1 : feat.End]
								if feat.Strand == '-' {
									fSeq = anti(fSeq)
								}

								tModelSc += model.Score(fSeq)
								p1 := len(fSeq) - model.K
								if p1 > 0 {
									tModelPos += p1
								}

								tBaseSc += baseModel.Score(fSeq)
								p2 := len(fSeq) - baseModel.K
								if p2 > 0 {
									tBasePos += p2
								}
							}
						}
					}
				}

				if tModelPos > 0 && tBasePos > 0 {
					diff := (tBaseSc / float64(tBasePos)) - (tModelSc / float64(tModelPos))
					if diff > 0.0 {
						badgenes = append(badgenes, BadGeneRes{Gid: gGroup.Gene.Fid, Tid: tx.TxFeat.Fid, Diff: diff})
					}
				}
			}
		}

		sort.Slice(badgenes, func(i, j int) bool { return badgenes[i].Diff > badgenes[j].Diff })
		f, _ := os.Create("badgenes.txt")
		defer f.Close()
		writer := bufio.NewWriter(f)
		fmt.Fprintln(writer, "GeneID\tTranscriptID\tScoreDiff")
		for _, bg := range badgenes {
			fmt.Fprintf(writer, "%s\t%s\t%.4f\n", bg.Gid, bg.Tid, bg.Diff)
		}
		writer.Flush()
		fmt.Printf("Analysis complete. Found %d anomalous transcripts. Results written to badgenes.txt\n", len(badgenes))

	case "auto":
		if len(cmdArgs) < 1 {
			fmt.Println("Usage: gene_qc auto <url> [-k <kmer_size>]")
			return
		}
		url := cmdArgs[0]
		k := 3
		if len(cmdArgs) == 3 && cmdArgs[1] == "-k" {
			k, _ = strconv.Atoi(cmdArgs[2])
		}

		fmt.Println("Initializing bulk pipeline...")
		jobs, err := crawlDatabase(url)
		if err != nil {
			fmt.Printf("Failed to crawl directory: %v\n", err)
			return
		}

		fmt.Printf("Found %d valid FASTA/GFF3 pairs to process.\n", len(jobs))

		var wg sync.WaitGroup
		// Buffered channel acts as a concurrency semaphore (max 4 concurrent jobs)
		sem := make(chan struct{}, 4)

		for _, job := range jobs {
			wg.Add(1)
			sem <- struct{}{} // Block here if there are already 4 items in the channel
			go processGenomeWorker(job, k, verbose, &wg, sem)
		}

		wg.Wait() // Wait for all workers to finish
		fmt.Println("Bulk pipeline complete! Check the 'bulk_out/' directory.")

	default:
		fmt.Printf("Unknown command: %s\n", command)
	}
}
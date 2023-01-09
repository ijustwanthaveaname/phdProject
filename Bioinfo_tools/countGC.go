package main

import (
	"bufio" // read lines from os.Open() one to one line
	"flag"  // read command
	"fmt"
	"io"      // io.EOF
	"os"      // open file
	"strings" // count gc
)

// func init() {
// 	flag.StringVar(&path, "path", "genome.fna", "path to genome file")
// }

func main() {
	var (
		path       string
		total      int
		g_count    int
		c_count    int
		n_count    int
		gc_content float64
	)
	flag.StringVar(&path, "path", "genome.fna", "path to genome file")
	flag.Parse()
	file, err := os.Open(path)
	if err != nil {
		fmt.Println(err)
	}
	defer file.Close()
	reader := bufio.NewReader(file)
	for {
		str, err := reader.ReadString('\n')
		if err == io.EOF {
			break
		}
		if str[0:1] != ">" {
			g_count += strings.Count(str, "G") + strings.Count(str, "g")
			c_count += strings.Count(str, "C") + strings.Count(str, "c")
			n_count += strings.Count(str, "N") + strings.Count(str, "n")
			total += len(str) - 1
		}
	}
	gc_content = float64(g_count+c_count) / float64(total-n_count)
	fmt.Printf("GC content is %.2f%% \n", gc_content*100)
}

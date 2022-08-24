###replace the same function in the fastq_filter.cpp
bool CFastqFilter::FilterRead()
{
        uint32 read_len = static_cast<uint32>(seq_desc.read_end - seq_desc.read_start);
        read.assign((char*)input_part + seq_desc.read_start, read_len);

        kmc_api.GetCountersForRead(read, counters);
        uint32 valid_kmers = 0;
        uint32 total_counter=0;
        for (auto counter : counters)
    {
        total_counter=total_counter+counter;
        printf("%d\n",counter);
        if (counter >= n_min_kmers && counter <= n_max_kmers)
        {
            ++valid_kmers;

        }
    }
        if (use_float_value)
        {
                uint32 min = static_cast<uint32>(f_min_kmers * (read_len - kmer_len + 1));
                uint32 max = static_cast<uint32>(f_max_kmers * (read_len - kmer_len + 1));
                if (valid_kmers >= min && valid_kmers <= max)
                        return true;
                return false;
        }
        else
        {
        float average =float(total_counter)/read_len;
        float percentage =float(valid_kmers)/read_len;
        printf("average is %.2f\n",average);
        printf("percentage is %.2f\n",percentage);
                if ( percentage > 0.9 )
        {
            printf("check true  %.2f\n",percentage);
                        return true;
        }
                return false;
        }
}

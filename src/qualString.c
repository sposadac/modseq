#include<pandaseq-plugin.h>
#include<stdlib.h>

HELP("Returns the original quality string", "qualString:filename");

VER_INFO("1.0");

struct data {
	PandaWriter writer;
};


static bool check_func(
	PandaLogProxy logger,
	const panda_result_seq *sequence,
	struct data *data) {
	size_t pos;
	int it = sequence->forward_length - sequence->overlap - 1;
	char ntFor;
	char ntRev;
	char ntSeq;
	char qFor;
	char qRev;

	panda_writer_append_c(data->writer, '@');
	panda_writer_append_id(data->writer, &sequence->name);
	panda_writer_append_c(data->writer, '\n');	
	for (pos = 0; pos < sequence->sequence_length; pos++) {
		panda_writer_append_c(data->writer, panda_nt_to_ascii(sequence->sequence[pos].nt));
	}
	panda_writer_append(data->writer, "\n+\n");
	for (pos = 0; pos < sequence->sequence_length; pos++) {
		qFor = (char) 0;
		qRev = (char) 0; 
		if(pos < sequence->forward_length) {
			ntFor = sequence->forward[pos].nt; 
			qFor = sequence->forward[pos].qual;
		}  
		else{ntFor = (char) 0;}
		if(it < pos || it < 0){
			if(sequence->reverse_length+it-pos >= 0) {
				ntRev = sequence->reverse[sequence->reverse_length+it-pos].nt;
				qRev = sequence->reverse[sequence->reverse_length+it-pos].qual;
			} 
		} else{ntRev = (char) 0;}
		ntSeq = sequence->sequence[pos].nt;		
		if(ntSeq == ntFor && ntSeq == ntRev){
			if(qFor >= qRev) {panda_writer_append_c(data->writer, 33 + qFor);} else if(qRev > qFor) {panda_writer_append_c(data->writer, 33 + qRev);}
		} 
		else if(ntSeq == ntFor) {
			panda_writer_append_c(data->writer, 33 + qFor);
		}
		else if(ntSeq == ntRev) {
			panda_writer_append_c(data->writer, 33 + qRev);
		}
	}
	panda_writer_append_c(data->writer, '\n');
	panda_writer_commit(data->writer);
	return true;
}

static void destroy_func(
	struct data *data) {
	panda_writer_unref(data->writer);
	free(data);
}

OPEN {
	struct data data;
	if (args == NULL || *args == '\0') {
		panda_log_proxy_write_str(logger, "Need a output file name.\n"); /*Print a string to the log.*/
		return false;
	}

	data.writer = panda_writer_open_file(args, false);
	if (data.writer == NULL) return false;
	
	
	*check = (PandaCheck) check_func;
	*precheck = NULL;
	*user_data = PANDA_STRUCT_DUP(&data);
	*destroy = (PandaDestroy) destroy_func; 
	return true;
}

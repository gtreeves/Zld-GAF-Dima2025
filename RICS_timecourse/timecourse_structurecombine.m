function data = timecourse_structurecombine(Data,NC_list)

v = find(~cellfun(@isempty,Data));
data = Data{v(1)}(1); % prealloc
xilimit = 64;
W_max = data.metadata.lsminf2.DimensionX;
fnames = fieldnames(data);
nnc = length(NC_list);
ntmax = 200; % comfortable max value
exceptions = {'pth','basename','genotype','side','modeltot',...
	'modelnuc','modelcyt','modelcc','modelnuc_ch'};

for i = 1:length(fnames)
	fnames1 = fnames{i};
	if ~any(strcmp(exceptions,fnames1))

		%
		% Structures (like "fittot") and char (like "filename")
		% will each be grouped into a cell variable across the NCs. The
		% nuclear and cytoplasmic logical masks also fit into this
		% category.
		%
		if isstruct(data.(fnames1)) ||  ischar(data.(fnames1)) || ...
				((isnumeric(data.(fnames1)) || islogical(data.(fnames1))) && ...
				size(data.(fnames1),2) > 1 && size(data.(fnames1),2) > W_max/2)
			c_1 = cell(1,nnc);
			for j = 1:nnc
				data1 = Data{j};
				for jj = 1:length(data1)
					c_1{j}{jj} = data1(jj).(fnames1);
				end
			end

		%
		% Numerical arrays that are either scalar or vector (like
		% "Anuc") and the datetime array (ie, "t") will each be grouped
		% into a ntmax-by-5 variable across the NCs.
		%
		elseif (isnumeric(data.(fnames1)) || isdatetime(data.(fnames1))) && ...
			(size(data.(fnames1),2) == 1 || size(data.(fnames1),1) == 1)
			if isnumeric(data.(fnames1))
				c_1 = NaN(ntmax,nnc);
			elseif isdatetime(data.(fnames1))
				c_1 = NaT(ntmax,nnc);
			end
			jjmax = 1;
			for j = 1:nnc
				data1 = Data{j};
				jj_idx = 1;
				for jj = 1:length(data1)
					c1 = data1(jj).(fnames1);
					c_1(jj_idx:(jj_idx+length(c1)-1),j) = c1;
					jj_idx = jj_idx+length(c1) + 1;
					jjmax = max(jjmax,jj_idx-2);
				end
			end
			c_1(jjmax+1:end,:) = [];



		%
		% Numerical arrays that are h-by-w-by-nt (ngroups actually), like
		% the ACFs, will be stored as numerical arrays, with the first two
		% dimensions being ntmax,nnc as normal, and (perhaps
		% unconventionally) the 3rd and 4th dimensions being h and w. We
		% are assuming we only need the first 64 elements of each ACF/CCF
		% in either direction (fast or slow direction). Note that the ACFs
		% and CCF are less or equal to one half the original image size
		% (ie, usually 1024-by-1024, or 512, or 256, or whatever...this
		% size is stored in W_max)
		%
		elseif isnumeric(data.(fnames1)) && size(data.(fnames1),2) > 1 && ...
				size(data.(fnames1),2) <= W_max/2
			
			c_1 = cell(1,nnc);

			for j = 1:nnc
				data1 = Data{j};
				for jj = 1:length(data1)
					
					c1 = data1(jj).(fnames1);
					c_1{j}{jj} = c1(1:xilimit,1:xilimit,:);
				end
			end

			

		end
		data.(fnames1) = c_1;

	end
end











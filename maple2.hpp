#pragma once

namespace libPixie
{
	using string = std::string;
	using strings = std::vector<std::string>;

	constexpr uint16 NOENTITY = 65535;

	struct Chunk
	{
		Float2 unit{ 4, 2 };			

		Float3 pos{ 0,0,0 };
		Float3 scale{ 1,1,1 };
		Float3 rotate{ 0,0,0 };

		Float3 motionpos{ 0,0,0 };
		Float3 motionscale{ 1,1,1 };
		Float3 motionrotate{ 0,0,0 };

		Mat4x4 mat()
		{
			Quaternion q1 = Quaternion::RollPitchYaw(ToRadians(motionrotate.x),
													 ToRadians(motionrotate.y),
													 ToRadians(motionrotate.z));
			Quaternion q2 = Quaternion::RollPitchYaw(ToRadians(rotate.x),
													 ToRadians(rotate.y),
													 ToRadians(rotate.z));
			Float3 sca = motionscale * scale;
			Float3 tra = motionpos + pos;

			Mat4x4 mat =
				Mat4x4::Identity().Scale(Float3{ -sca.x,sca.y,sca.z }) * Mat4x4(q1 * q2) *
				Mat4x4::Identity().Translate(tra);

			return mat;
		}
	};

	struct Entity
	{
		std::string meshtag;	 
		int32 meshid;			

		Float3 pos{ 0,0,0 };	
		Float3 scale{ 1,1,1 };	
		Quaternion qrotate = Quaternion::Identity();
	};

	bool ascZ(const Entity& L, const Entity& R)
	{
		return (L.pos.z == R.pos.z) ? (L.pos.x < R.pos.x) : (L.pos.z < R.pos.z);
	}

	bool ascX(const Entity& L, const Entity& R)
	{
		return (L.pos.x == R.pos.x);
	}

	struct Maple
	{
		std::string grid = {};		 
		Size size = {};				 
		Array<uint16> refgrid;		 
		Array<Entity> displaylist;	 

		Maple() = default;
		Maple(const std::string& map) : grid{ map }
		{
			size = GetSize(grid);	
			grid.erase(std::remove(grid.begin(), grid.end(), '\n'), grid.end());		
		}

		Size GetSize(std::string& grid)
		{
			int w = grid.find('\n');
			int h = std::count(grid.cbegin(), grid.cend(), '\n');
			return Size{ w, h };
		}

		static size_t GetPoint(Point& point, size_t start, std::string& mapstr, std::string& tag, uint32 width)
		{
			size_t index = mapstr.find(tag, start);
			if (index == std::string::npos) return std::string::npos;	
			point.x = index % width;
			point.y = index / width;
			return index;
		}

		char GetChar(const Float2& pos, const OrientedBox& ob) const
		{
			Float3 gridpos = {};
			return GetChar(pos, ob, gridpos);
		}
		char GetChar(const Float2& newpos, const OrientedBox& ob, Float3& gridpos) const
		{
			if (size == Size{ 0,0 }) assert(1);
			Float2 unit{ newpos.x / ob.w , newpos.y / ob.d };

			gridpos.x = (int)Clamp((float)size.x * unit.x, 0.0f, (float)size.x - 1);
			gridpos.y = (int)Clamp((float)size.y * unit.y, 0.0f, (float)size.y - 1);
			gridpos.z = gridpos.y * size.x + gridpos.x;

			return grid.at(gridpos.z);
		}

		void Replace(std::string& text, const std::string& from, const std::string& to) const
		{
			if (!text.empty())
			{
				std::string::size_type pos = 0;
				while ((pos = text.find(from, pos)) != std::string::npos)
				{
					text.replace(pos, from.length(), to);
					pos += to.length();
				}
			}
		}

	};

	Mat4x4 MatTSR(Float3 tra = { 0,0,0 }, Float3 scale = { 1,1,1 }, Quaternion qrot = Quaternion::Identity())
	{
		return Mat4x4::Identity().Scale(Float3{ -scale.x,scale.y,scale.z }) *
			Mat4x4(qrot) * Mat4x4::Identity().Translate(tra);
	}

	void initMaple(Maple& maple, Array<Entity>& entities, Size unit = Size{ 4,2 })
	{
		if (maple.grid.size()) assert(1);

		if (maple.displaylist.size() == 0)
		{
			int CENTER = maple.size.x / (4 * 2);
			for (auto& ent : entities)	
			{
				size_t index = 0;
				for (;;)
				{
					Point point = {};
					index = Maple::GetPoint(point, index, maple.grid, ent.meshtag, maple.size.x);
					if (index == std::string::npos) break;
					ent.pos = Float3{ point.x / unit.x, 0, point.y / unit.y };	
					maple.displaylist.emplace_back(ent);		
					index++;
				}
			}

			auto& DL = maple.displaylist;
			std::sort(DL.begin(), DL.end(), ascZ);

			maple.size.x /= unit.x;
			maple.size.y /= unit.y;
			maple.refgrid.resize(maple.size.x * maple.size.y, NOENTITY);

			for (uint i = 0; i < DL.size(); i++)
			{
				maple.refgrid[(DL[i].pos.z * maple.size.x) + DL[i].pos.x] = i;
				DL[i].pos.x -= CENTER;
			}
		}
	}
	void drawMaple(PixieMesh& mesh, Maple& maple, Mat4x4 mat,
					double progress, LineString3D curve, Rect process = Rect{ 0,0,-1,0 })
	{
		if (maple.displaylist.size() == 0) assert(1);

		Float3 sca, pos;
		Quaternion qrot;
		mat.decompose(sca, qrot, pos);

		if (process.w == -1) process.size = maple.size;
		else
		{
			Clamp(process.x, 0, maple.size.x);
			Clamp(process.y, 0, maple.size.y);
			Clamp(process.w, 0, maple.size.x);
			Clamp(process.h, 0, maple.size.y);
		}

		int ww = maple.size.x;
		process.y = curve.size() * progress;
		for (int y = 0; y < process.h; y++)
		{
			int yy = process.y + y;
			yy = (yy > curve.size()) ? (yy - curve.size() - 1) : yy;

			for (int x = 0; x < process.w; x++)
			{
				int xx = process.x + x;
				xx = (xx > maple.size.x) ? (xx - maple.size.x - 1) : xx;

				uint16 ref = maple.refgrid[yy * ww + xx];
				if (ref < maple.displaylist.size())
				{
					Entity& DL = maple.displaylist[ref];
					Mat4x4 mat = MatTSR(DL.pos + pos, DL.scale * sca, DL.qrotate * qrot);
					mesh.drawMesh(mat, DL.meshid);
				}
			}
		}
	}

}

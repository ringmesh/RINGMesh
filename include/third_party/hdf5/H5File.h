// C++ informative line for the emacs editor: -*- C++ -*-
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef __H5File_H
#define __H5File_H

#ifndef H5_NO_NAMESPACE
namespace H5 {
#endif

/*! \class H5File
    \brief Class H5File represents an HDF5 file.

    It inherits from H5Location and CommonFG.
*/
class H5_DLLCPP H5File : public H5Location, public CommonFG {
   public:
	// Creates or opens an HDF5 file.
	H5File( const char* name, unsigned int flags,
	   const FileCreatPropList& create_plist = FileCreatPropList::DEFAULT,
	   const FileAccPropList& access_plist = FileAccPropList::DEFAULT );
	H5File( const H5std_string& name, unsigned int flags,
	   const FileCreatPropList& create_plist = FileCreatPropList::DEFAULT,
	   const FileAccPropList& access_plist = FileAccPropList::DEFAULT );

	// Open the file
	void openFile(const H5std_string& name, unsigned int flags,
	    const FileAccPropList& access_plist = FileAccPropList::DEFAULT);
	void openFile(const char* name, unsigned int flags,
	    const FileAccPropList& access_plist = FileAccPropList::DEFAULT);

	// Close this file.
	virtual void close();

	// Gets the access property list of this file.
	FileAccPropList getAccessPlist() const;

	// Gets the creation property list of this file.
	FileCreatPropList getCreatePlist() const;

	// Retrieves the file size of an opened file.
	hsize_t getFileSize() const;

	// Returns the amount of free space in the file.
	hssize_t getFreeSpace() const;

	// Returns the number of opened object IDs (files, datasets, groups
	// and datatypes) in the same file.
	ssize_t getObjCount(unsigned types = H5F_OBJ_ALL) const;

	// Retrieves a list of opened object IDs (files, datasets, groups
	// and datatypes) in the same file.
	void getObjIDs(unsigned types, size_t max_objs, hid_t *oid_list) const;

	// Returns the pointer to the file handle of the low-level file driver.
	void getVFDHandle(void **file_handle) const;
	void getVFDHandle(const FileAccPropList& fapl, void **file_handle) const;
	void getVFDHandle(FileAccPropList& fapl, void **file_handle) const; // kept for backward compatibility

	// Determines if a file, specified by its name, is in HDF5 format
	static bool isHdf5(const char* name );
	static bool isHdf5(const H5std_string& name );

	// Reopens this file.
	void reOpen();	// added for better name

#ifndef DOXYGEN_SHOULD_SKIP_THIS
	void reopen();  // obsolete in favor of reOpen()

	// Gets the file id
	virtual hid_t getLocId() const;

	// Creates an H5File using an existing file id.
	H5File(hid_t existing_id);

#endif // DOXYGEN_SHOULD_SKIP_THIS

	///\brief Returns this class name.
	virtual H5std_string fromClass () const { return("H5File"); }

	// Throw file exception.
	virtual void throwException(const H5std_string& func_name, const H5std_string& msg) const;

	// Default constructor
	H5File();

	// Copy constructor: makes a copy of the original H5File object.
	H5File(const H5File& original);

	// Gets the HDF5 file id.
	virtual hid_t getId() const;

	// H5File destructor.
	virtual ~H5File();

   protected:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
	// Sets the HDF5 file id.
	virtual void p_setId(const hid_t new_id);
#endif // DOXYGEN_SHOULD_SKIP_THIS

   private:
	hid_t id;	// HDF5 file id

	// This function is private and contains common code between the
	// constructors taking a string or a char*
	void p_get_file( const char* name, unsigned int flags, const FileCreatPropList& create_plist, const FileAccPropList& access_plist );

};
#ifndef H5_NO_NAMESPACE
}
#endif
#endif // __H5File_H
